#!/usr/bin/env perl
# Sam Shepard - 2019.03

use strict;
use warnings;

use Getopt::Long;
my ($message,$inPlace,$strandField);
GetOptions('in-place|I' => \$inPlace, 'strand-field|S' => \$strandField );

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <fasta> <IDs> [options]\n";
	$message .= "\t\t-I|--in-place\t\tWrite out file in place rather than to STDOUT.\n";
	$message .= "\t\t-S|--strand-field\tExpects tab-delimited fields (ID:first, strand:last). If strand is '-', reverse complement.\n";
	die($message."\n");
}

$/ = "\n";
my %IDs = (); 
my $line = '';
my @f = ();
open(TAB,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while($line = <TAB>) {
	chomp($line);
	@f = split("\t",$line);
	if ( defined($f[0]) && $f[0] ne '' ) {
		if ( !defined($strandField) || $f[$#f] eq '-' ) {
			$IDs{$f[0]} = 1;
		}
	}
}
close(TAB);

$/ = '>'; 
my ($seqID,$sequence,$record) = ('','','');
my $length;
my @lines = ();
if ( defined($inPlace) ) {
	open(FASTA,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
	my @records = <FASTA>; chomp(@records);
	close(FASTA);
	open(FASTA,'>',$ARGV[0]) or die("Cannot open $ARGV[0] for writing.\n");
	foreach $record ( @records ) {
		@lines = split(/\r\n|\n|\r/, $record);
		$seqID = shift(@lines);
		$sequence = lc(join('',@lines));

		$length = length($sequence);
		if ( $length == 0 ) { next; }

		if ( defined($IDs{$seqID}) ) {
			$sequence = reverse( $sequence );
			$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
		}

		print FASTA '>',$seqID,"\n",$sequence,"\n";
	}
} else {
	open(FASTA,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
	while( $record = <FASTA> ) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$seqID = shift(@lines);
		$sequence = lc(join('',@lines));

		$length = length($sequence);
		if ( $length == 0 ) {
			next;	
		}

		if ( defined($IDs{$seqID}) ) {
			$sequence = reverse( $sequence );
			$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
		}

		print STDOUT '>',$seqID,"\n",$sequence,"\n";
	}
	close(FASTA);
}
