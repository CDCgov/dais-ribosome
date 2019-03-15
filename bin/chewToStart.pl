#!/usr/bin/env perl

use Getopt::Long;
GetOptions(	'ref|R=s' => \$referenceFile, 'length|L=i' => \$referenceLength, 'flex-length|X=i' => \$flexLen
		);


if ( -t STDIN && scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <nts.fasta> [options]\n";
	$message .= "\t\t-R|--ref\t\tReference file name.\n";
	$message .= "\t\t-L|--length\t\tReference length.\n";
	$message .= "\t\t-X|--flex-length\tAmount the 'chewed' sequence may be less than the reference.\n";
	die($message."\n");
}

$/ = '>';
if ( defined($referenceFile) ) {
	open(REF,'<',$referenceFile) or die("Cannot open $referenceFile for reading.\n");
	while ( $record = <REF> ) {
		chomp($record);
		@lines = split(/\r\n|\n|\r/, $record);
		$id = shift(@lines);
		$sequence = uc(join('',@lines));
		$length = length($sequence);
		if ( $length == 0 ) {
			next;
		} else {
			$REF_LEN = $length;
			last;
		}
	}
} elsif ( defined($referenceLength) && $referenceLength > 0 ) {
	$REF_LEN = $referenceLength;
} else {
	$REF_LEN = 0;
}

if ( defined($flexLen) ) { $REF_LEN -= $flexLen; }
while ( $record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$id = shift(@lines);
	$sequence = uc(join('',@lines));
	$length = length($sequence);
	if ( $length == 0 ) { next; }

	if ( $REF_LEN < $length && $sequence =~ /ATG/ ) {
		$newSeq = substr($sequence,$-[0]);
		$newLen = length($newSeq);
		if ( $newLen >= $REF_LEN ) {
			print '>',$id,"\n",$newSeq,"\n";
		} else {
			print '>',$id,"\n",$sequence,"\n";	
		}
	} else {
		print '>',$id,"\n",$sequence,"\n";	
	}
}
