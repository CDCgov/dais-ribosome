#!/usr/bin/env perl

if ( scalar(@ARGV) != 2 ) {
	$message = "perl $0 <config_file> <gene_segment_subtype>\n";
	die($message."\n");
}

$filename= $ARGV[0];
$gene_selected = $ARGV[1];

$/ = "\n";
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
$default = <IN>; chomp($line); $gene = '';
while($line=<IN>) {
	chomp($line);
	($gene,$params,$ref) = split("\t",$line);
	if ( index($gene_selected, $gene) > -1 ) {
		($match,$mismatch,$gapopen,$gapextend) = split(' ',$params);
		if ( $match < 1 ) { $match = 1; }
		if ( $mismatch < 1 ) { $mismatch = 1; }
		if ( $gapopen < 1 ) { $gapopen = 1; }
		if ( $gapextend < 1 ) { $gapextend = 1; }

		print " -m $match -x $mismatch -o $gapopen -e $gapextend ";
		exit;
	}
}


($gene,$params,$ref) = split("\t",$default);
($match,$mismatch,$gapopen,$gapextend) = split(' ',$params);
if ( $match < 1 ) { $match = 1; }
if ( $mismatch < 1 ) { $mismatch = 1; }
if ( $gapopen < 1 ) { $gapopen = 1; }
if ( $gapextend < 1 ) { $gapextend = 1; }

print " -m $match -x $mismatch -o $gapopen -e $gapextend ";

