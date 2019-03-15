#!/usr/bin/env perl
# Sam Shepard

use File::Basename;
use Getopt::Long;
GetOptions(	'gene-segment|G=s' => \$geneSegment, 'reference-id|R=s' => \$referenceID
		);

if ( scalar(@ARGV) != 3 ) {
	die("Usage:\n\t$0 <fasta> <product-table> <prefix> [-G|--gene-segment <CTYPE>] [-R|--reference-id <STR>\n");
}

$PROG = basename($0,'.pl');

$/ = "\n"; $max = 0; $productsFound = 0;
open(PROD,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while($line=<PROD>) {
	chomp($line);
	($segment,$peptide,$headerInfo,$coords,$prefix,$suffix) = split("\t",$line);
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
	$context{$peptide} = [lc($prefix),lc($suffix)];
	for( $i=0;$i<scalar(@coordList);$i++ ) {
		($start,$stop) = split(',',$coordList[$i]);
		$index = $start - 1;
		$L = $stop - $index;
		if ( $stop > $max ) {
			$max = $stop;
		}
		$exons{$peptide}[$i][0] = $index;
		$exons{$peptide}[$i][1] = $L;
	}
}
close(PROD);

if ( defined($geneSegment) && $productsFound == 0 ) {
	die("$PROG:\tNo products found for $geneSegment.\n");
}

%handles = ();
@peptides = sort { $a cmp $b } keys(%fields);
#foreach $peptide ( @peptides ) {
#	$filename = $ARGV[2] . '-' . $peptide . '.fasta';
#	open($peptide,'>', $filename) or die("Cannot open $filename for writing.\n");
#}

$filename = $ARGV[2] .'.products';
open(OUT,'>', $filename) or die("Cannot open $filename for writing.\n");
$/ = '>';
open(FASTA,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while ( $record = <FASTA> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$id = shift(@lines);
	$sequence = lc(join('',@lines));
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	} elsif ( $max > $length ) {
		die("Found a sequence shorter ($length) than the last coordinate position ($max).\n");
	} else {
		foreach $p ( @peptides ) {
			($prefix,$suffix) = @{$context{$p}};
			$cds = '';
			for($i=0;$i<scalar(@{$exons{$p}});$i++) {
				$cds .= substr($sequence,$exons{$p}[$i][0],$exons{$p}[$i][1]);
			}

			if ( $prefix ne "" && $prefix ne substr($sequence,0,length($prefix)) ) {
				next;
			} elsif ( $suffix ne "" && $suffix ne substr($sequence, -length($suffix)) ) {
				next;
			}

			$length = length($cds);
			if ( $length % 3 != 0 ) { die("$PROG:\tNot in triplets ($length) for peptide '$p'.\n"); }
			print OUT '>',$id,'|',$fields{$p},"\n",$cds,"\n";
		}
	}
}
close(FASTA);
close(OUT);
#foreach $peptide ( @peptides ) {
#	close($handles{$peptide});
#}
