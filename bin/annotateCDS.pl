#!/usr/bin/perl
# Samuel Shepard - 6.2018
# Annotate the CDS of sequences in DAIS

use File::Basename;

if (  scalar(@ARGV) < 5 ) {
	$message = "\nUsage:\n\t$0 <0-spec> <2-original-fasta> <5-product-ins-c> <5-aligned-fas-c> <4-aligned-fas/ins ...> \n";
	die($message."\n");
}

$specfile=shift(@ARGV);
$origfile=shift(@ARGV);
$insprodfile=shift(@ARGV);
$fasprodfile=shift(@ARGV);

# FUNCTIONS #
sub removeElongation($) {
	my $seq = $_[0];
	$seq =~ s/[A-Z]+$//;
	return $seq;
}

sub sequenceByCigar($$$) {
	my ($sequence,$cigar,$offset) = ($_[0],$_[1],$_[2]);
	my ($result,$state) = ('','');
	my $length = 0;
	while( $cigar =~ m/(\d+)([MIDSHN])/g ) {
		($length,$state) = ($1,$2);
		if ( $state =~ /[MI]/ ) {
			$result .= substr($sequence,$offset,$length);
			$offset += $length;
		}
	}
	return $result;
}

sub getRef($) {
	my $file = basename($_[0]);
	my @pieces = split('\.',$file);
	return $pieces[0];
}


sub getSubstringOffset($$) {
	my ($original, $alignable) = (lc($_[0]),lc($_[1]));
	my $leftpad = 0;
	if ( $alignable =~ /^(\.+)/ ) {
		$leftpad = length($1);		
	}
	$alignable =~ tr/.-//d;

	if ( $original =~ m/\Q$alignable\E/ ) {
		return ($-[0] - $leftpad);
	} else {
		print STDERR "ERROR, alignable not found.\n>Original\n$original\n>Alignable\n$alignable\n\n";
		return 0.5;
	}
}

sub getSubstringCoords($$) {
	my ($original, $alignable) = (lc($_[0]),lc($_[1]));
	$alignable =~ tr/.-//d;

	if ( $original =~ m/\Q$alignable\E/ ) {
		return ($-[0],$+[0]);
	} else {
		print STDERR "ERROR, alignable not found.\n>Original\n$original\n>Alignable\n$alignable\n\n";
		return (-1,-1);
	}
}

sub sequenceToCigar($) {
	my $seq = $_[0];
	$seq =~ tr/A-Z/I/;
	$seq =~ tr/a-z/M/;
	$seq =~ tr/-/D/;
	$seq =~ tr/./N/;

	my $cigar = '';
	while( $seq =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/g ) {
		$cigar .= length($1).substr($1,0,1);
	}
	return $cigar;
}

sub condenseCigar($) {
	my $cig = $_[0];
	my $cigar = '';
	my $state = '';
	while( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/g ) {
		$state = $1;
		$cigar .= length($state);
		$cigar .= substr($state,0,1);
	}
	return $cigar;
}

sub sequenceToStates($) {
	my $seq = $_[0];
	$seq =~ tr/A-Z/I/;
	$seq =~ tr/a-z/M/;
	$seq =~ tr/-/D/;
	$seq =~ tr/./N/;
	return $seq;
}

sub addInsertions($$) {
	my $seq = lc($_[0]);
	my $inserts = $_[1];
	my ($offset,$pos) = 0;
	my $insert = '';
	foreach $pos ( sort { $a <=> $b } keys(%{$inserts}) ) {
		$insert = $inserts->{$pos};
		substr($seq,int($pos)+$offset,0) = uc($insert);
		$offset += length($insert);
	}

	return $seq;
}

sub addInsertionsBounded($$$) {
	my $seq = lc($_[0]);
	my $inserts = $_[1];
	my $offset = $_[2];
	my $pos = 0;
	my $L = length($seq);
	foreach $pos ( sort { $a <=> $b } keys(%{$inserts}) ) {
		# 1 - based check
		if ( ($pos+$offset) > $L || ($pos+$offset) < 1 ) { last; }
		substr($seq,$pos+$offset,0) = uc($inserts->{$pos});
		$offset += length($inserts->{$pos});
		$L = length($seq);
	}

	return $seq;
}

# process parameters
$PROG = basename($0,'.pl');

# process specifications
$/ = "\n"; %specs = (); %segmentByRefPep = ();
$max = 0; $productsFound = 0;
open( SPEC, '<', $specfile ) or die("$PROG ERROR: Could not open $specfile for reading.\n");
while($line=<SPEC>) {
	chomp($line);
	($segment,$peptide,$headerInfo,$coords,$prefix,$suffix) = split("\t",$line);
	($ref_id,$peptide2) = split('\|',$headerInfo);
	@coordList = split(';',$coords);

	$segmentByRefPep{$ref_id}{$peptide} = $segment;
	for( $i=0;$i<scalar(@coordList);$i++ ) {
		($start,$stop) = split(',',$coordList[$i]);
		$specs{$ref_id}{$peptide}[$i][0] = $start - 1;
		$specs{$ref_id}{$peptide}[$i][1] = $stop - $start + 1;
	}
}
close(SPEC);

# process original fasta
open( ORIG, '<', $origfile ) or die("$PROG ERROR: Could not open $origfile for reading.\n");
$/ = ">"; %originals = ();
while( $record = <ORIG> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$seq = lc(join('',@lines));

	$length = length($seq);
	if ( $length == 0 ) { 
		next;
	} else {
		@id = split('\|',$header);
		$originals{$id[0]} = $seq;
	}
}
close(ORIG);

# process alignment (gene segment level) insertions
$/ = "\n"; %segmentInsertions = ();
foreach $file (@ARGV) {
	if ( $file =~ /\.ins\.txt$/ ) {
		$ref_id = getRef($file);
	} else {
		next;
	}

	open(INS, '<', $file) or die("$PROG ERROR: could not open $file for reading.\n");
	while($line = <INS> ) {
		chomp($line);
		($compound_id,$upstream_position,$insert) = split("\t",$line);
		($flu_seq_id,$segment) = split('\|',$compound_id);
		$segmentInsertions{$ref_id}{$segment}{$flu_seq_id}{$upstream_position} = $insert;
	}
	close(INS);
}


# [INPUT]: "240934|A_HA_H3|HK4801|PB1	2274	TGA"
# process product (nucleotide gene transcript) insertions
$/ = "\n"; %productInsertions = ();
open(INS, '<', $insprodfile) or die("$PROG ERROR: could not open $insprodfile for reading.\n");
while($line = <INS> ) {
	chomp($line);
	($compound_id, $upstream_position, $insert) = split("\t",$line);
	($flu_seq_id,$segment,$ref_id,$peptide) = split('\|',$compound_id);
	$productInsertions{$ref_id}{$peptide}{$flu_seq_id}{$upstream_position} = $insert;
}
close(INS);


# process the segment coordinates to create a bounding box for segment alignments via reference
$/ = ">"; %segmentOffset = ();
foreach $file ( @ARGV ) {
	if ( $file =~ /fasta$/ ) {
		$ref_id = getRef($file);
	} else {
		next;
	}

	$total = $found = 0;
	open( ALIGNED, '<', $file ) or die("$PROG ERROR: Could not open $file for reading.\n");
	while( $record = <ALIGNED> ) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$header = shift(@lines);
		$seq = lc(join('',@lines));
		($flu_seq_id,$segment) = split('\|',$header);

		$length = length($seq);
		if ( $length == 0 ) { 
			next;
		} elsif ( defined($originals{$flu_seq_id}) ) {
			$found++;

			# Add back in insertions so the reference-aligned sequence will map to the original.
			# Remove the trailing insertions (3' elongation) if applicable.
			# Goal is to get the Reference boundary coords, which excludes elongation.
			if ( defined($segmentInsertions{$ref_id}{$segment}{$flu_seq_id}) ) {
				$seq = removeElongation(addInsertions($seq, \%{$segmentInsertions{$ref_id}{$segment}{$flu_seq_id}}))
			}	

			# Reference-Aligned Query to Original	
			$offset = getSubstringOffset($originals{$flu_seq_id},$seq);
			if ( $offset == 0.5 ) { 
				die("Issue with $segment / $ref_id / #$flu_seq_id!\n");
			}
			$segmentOffset{$segment}{$ref_id}{$flu_seq_id} = $offset;
		} else {
			print STDERR "Original pair not found: $flu_seq_id ( $segment / $ref_id )\n";
		}

		$total++;
	}
	close( ALIGNED );
}

# Sample header: >251324|A_NA_N2|HK4801|NA
# process final products and create a map between reference coordinates (codon number) and original nucleotide coordinates
$/ = ">";
open( PRODUCTS, '<', $fasprodfile ) or die("$PROG ERROR: Could not open $fasprodfile for reading.\n");
while( $record = <PRODUCTS> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$compound_id = shift(@lines);
	$seq_prod = lc(join('',@lines)); $seq_len = length($seq_prod);
	($flu_seq_id,$segment,$ref_id,$peptide) = split('\|',$compound_id);
	#$segment = $segmentByRefPep{$ref_id}{$peptide};

	$length = length($seq_prod);
	if ( $length == 0 ) { 
		next;
	} elsif ( defined($originals{$flu_seq_id}) ) {
		@exons = @{$specs{$ref_id}{$peptide}};
		$original = $originals{$flu_seq_id};

		# Aligned segment to original sequence offset. (1-based).
		$oriOffset = $segmentOffset{$segment}{$ref_id}{$flu_seq_id};
		($oriCoords, $pepCoords) = ('','');

		$pepOffset = 0; $first = 1;
		for($i = 0;$i < scalar(@exons);$i++ ) {
			($idx,$L) = @{$exons[$i]};

			# Get the cigar for the current exon, do not add insertions outside the valid bounds.
			# We process the exons from the peptide / product sequence.
			# I use the peptide term interchangably with product, although the product files are untranslated in DAIS terms.
			$exon = substr($seq_prod,$pepOffset,$L);
			$exonCigar = sequenceToCigar( addInsertionsBounded($exon, \%{$productInsertions{$ref_id}{$peptide}{$flu_seq_id}},-$pepOffset) );

			# Let $idx + 1 = offset for the peptide within the segment alignment
			# Thus the starting original coordinate adds the peptide to ref and ref to original offsets. 
			$oriCursor = $idx + 1 + $oriOffset;

			# Peptide cursor is after previous exon lengths (offsets). Start at 1.
			$pepCursor = $pepOffset+1;		

			# Create a map from the original coordinate space to the peptide product coordinate space.
			# Based on the OPeration, we advance each coordinate system relative to the other.
			while ( $exonCigar =~ m/(\d+)([MDNI])/g ) {
				($inc,$op) = ($1,$2);
				if ( $op eq 'N' ) {
					if  ( $first ) {
					#	$pepCoords .= $pepCursor.'..';
						$pepCursor += $inc;
					#	$pepCoords .= ($pepCursor-1).';';
					#	$oriCoords .= ($oriCursor-1).';';
						$first = 0;
					}
				} elsif ( $op eq 'M' ) {
					$pepCoords .= $pepCursor.'..';
					$oriCoords .= $oriCursor.'..';
					
					$pepCursor += $inc;
					$oriCursor += $inc;

					$pepCoords .= ($pepCursor-1).';';
					$oriCoords .= ($oriCursor-1).';';
				} elsif ( $op eq 'I' ) {
					$pepCoords .= ($pepCursor-1).';';
					$oriCoords .= $oriCursor.'..';

					$oriOffset += $inc;	# Take into account insertions in the original sequence.
					$oriCursor += $inc;

					$oriCoords .= ($oriCursor-1).';';
				} elsif ( $op eq 'D' ) {
#					$oriCoords .= ($oriCursor-1).';';
#					$pepCoords .= $pepCursor.'..';

					$pepCursor += $inc;
					$oriOffset -= $inc;	# Likewise deletions are removed from the offset.

#					$pepCoords .= ($pepCursor-1).';';
				} else {
					die("$op : Unknown\n");
				}
			}
			$pepOffset += $L;	# reflect that the exon has moved forward
		}

		chop($pepCoords);chop($oriCoords);
		print join('|',($flu_seq_id,$segment,$ref_id,$peptide)),"\t",$oriCoords,"\t",$pepCoords,"\n";
	} else {
		print STDERR "Original pair not found: $flu_seq_id ( $segment / $peptide / $ref_id )\n";
	}
}
close(PRODUCTS);
exit(0);
