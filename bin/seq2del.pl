#!/usr/bin/env perl
# Sam Shepard - 2020-06

use POSIX;
use strict;
use warnings;

use Getopt::Long;
my $genomeDeletions = 0;
GetOptions( 'genome|G' => \$genomeDeletions
	);

use constant { 
	AA_ALN 	=> 6,
	CDS_ID	=> 7,
	CDS_ALN	=> 11,
	GEN_ALN	=> 7
};


my $aln = CDS_ALN;
my $pre_end = 4;
if ( $genomeDeletions ) {
	$aln = GEN_ALN;
	$pre_end = 3;
}

if ( scalar(@ARGV) != 1 ) {
	die("\nUsage:\n\t$0 <input.seq> [--genome|-G]\n\n");
}

$/ = "\n";
open(IN,'<',$ARGV[0]) or die("$0 ERROR: cannot open $ARGV[0] for reading!\n");
while ( my $line = <IN> ) {
	chomp($line);
	my @f = split("\t",$line);

	if ( !defined($f[$aln]) ) {
		next;
	}
	
	my $prefix = join("\t",@f[0 .. $pre_end]);
	while( $f[$aln] =~ m/(-+)/g ) {
		my ($b, $e, $l) = ($-[0]+1, $+[0], $+[0]-$-[0]);
		my $in_frame = $b % 3 == 1 && $e % 3 == 0 ? 'true' : 'false';
		my ($aa_b, $aa_e, $aa_l ) = ( int( ($b-1)/3 + 1), int(($e-1) / 3 + 1), ceil($l/3) );


		if ( $genomeDeletions ) {
			print STDOUT $prefix,"\t",$b,"\t",$e,"\t",$l,"\n";
		} else {
			print STDOUT $prefix,"\t",$aa_b,"\t",$aa_e,"\t",$aa_l,"\t",$in_frame;
			print STDOUT "\t",$f[CDS_ID],"\t",$b,"\t",$e,"\t",$l,"\n";
		}
	}
}
