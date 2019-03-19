#!/usr/bin/env perl
# Sam Shepard

use File::Basename;
use Getopt::Long;
GetOptions(	'gene-segment|G=s' => \$geneSegment, 'reference-id|R=s' => \$referenceID
	);

if ( scalar(@ARGV) != 3 ) {
	die("Usage:\n\t$0 <ins-table> <product-table> <prefix> [-G|--gene-segment]\n");
}

$PROG = basename($0,'.pl');

$insertionTable = $ARGV[0];
$productTable = $ARGV[1];
$prefix = $ARGV[2];

$/ = "\n"; $max = $productsFound = 0;
%exons = %fields = %pMax = (); 
open(PROD,'<',$productTable) or die("Cannot open $productTable for reading.\n");
while($line=<PROD>) {
	chomp($line);
	($segment,$peptide,$headerInfo,$coords) = split("\t",$line);
	if ( $headerInfo =~ /\|/ ) {
		($refID,$junk) = split(/\|/,$headerInfo);
	} else {
		$refID = $headerInfo;
	}

	if ( defined($geneSegment) && $segment ne $geneSegment ) {
		next;
	} elsif ( defined($referenceID) && $refID ne $referenceID ) {
		next;
	} else {
		$productsFound++;
	}

	@coordList = split(';',$coords);
	$fields{$peptide} = $headerInfo;
	for( $i=0;$i<scalar(@coordList);$i++ ) {
		($start,$stop) = split(',',$coordList[$i]);

		if ( $max < $stop ) { $max = $stop; }
		if ( $pMax{$peptide} < $stop ) { $pMax{$peptide} = $stop; }
		$exons{$peptide}[$i][0] = $start;
		$exons{$peptide}[$i][1] = $stop;
	}
}
close(PROD);

if ( defined($geneSegment) && $productsFound == 0 ) {
	die("$PROG:\tNo products found for $geneSegment.\n");
}

@peptides = sort { $a cmp $b } keys(%fields);

$filename = $prefix.'.ins';
open(OUT,'>', $filename) or die("Cannot open $filename for writing.\n");
$/="\n";
open(INS,'<',$insertionTable) or die("Cannot open $insertionTable for reading.\n");
@lines = <INS>; chomp(@lines);
foreach $line ( @lines ) {
	@fields = split("\t",$line);
	if ( scalar(@fields) != 3 ) {
		die("Expected 3 fields in this format:\n\tID<tab>POSITION<tab>INSERT\n");
	} else {
		($id,$pos,$insert) = @fields;
		if ( $pos > $max ) {
			next;
			die("Invalid insertion detected: $pos >= $max\n");
		}

		foreach $p ( @peptides ) {
			$offset = 0; $newPos = 0;
			for($i=0;$i<scalar(@{$exons{$p}});$i++) {
				($start,$stop) = @{$exons{$p}[$i]};
				if ( $start <= $pos && $pos < $stop ) {
					$newPos = $pos - $start + 1 + $offset;
					print OUT $id,'|',$p,"\t",$newPos,"\t",$insert,"\n";
					$offset += ($stop-$start+1);
					last;
				}
				$offset += ($stop-$start+1);
			}
			if ( $pos == $max && $pos == $pMax{$p} ) {
				$newPos = $offset;
				print OUT $id,'|',$p,"\t",$newPos,"\t",$insert,"\n";
			}
		}
	}
}
close(INS);
close(OUT);
