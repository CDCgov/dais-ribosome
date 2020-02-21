#!/usr/bin/env perl
# Sam Shepard
# spec2db
# For outputing the DAIS ribosome spec/refs to a hadoop-ready table

use Digest::SHA qw(sha1_hex);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <spec> <refs>\n";
	die($message."\n");
}

sub nt_id2($) {
	my $seq = defined($_[0]) ? uc($_[0]) : '';
	$seq =~ tr/ :.~-//d;
	if ( $seq ne '' ) {
		return (sha1_hex($seq),$seq);
	} else {
		return ('\N','\N');
	}
}

my %seqByRefSeg = ();
$/ = '>'; 
open(FASTA,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while($record = <FASTA>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = uc(join('',@lines));

	($refID,$seg,$id) = (split('\|',$header))[0..2];
	if ( $id =~ /[A-Za-z]/ || length($sequence) == 0 ) {
		next;
	}

	if ( $refID ne '' ) {
		$seqByRefSeg{$refID}{$seg} = $sequence;
	}
}
close(FASTA);

$/ = "\n"; $extra = '';
open(SPEC,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($line = <SPEC> ) {
	chomp($line);
	($ctype,$prot,$refID,$prot2,$range_list) = split("\t|[|]",$line);
	$L = 0;
	@ranges = split(';',$range_list);
	foreach my $pair ( @ranges ) {
		($from,$to) = split(',',$pair);
		$L += $to - $from + 1;
	}
	$range_list =~ s/,/../g;
	$range_list =~ tr/;/,/;

	if ( length($seqByRefSeg{$refID}{$ctype}) > 0) {
		$extra = "\t" . join("\t",nt_id2($seqByRefSeg{$refID}{$ctype}));
	} else {
		$extra = '';
	}
	print STDOUT $ctype,"\t",$refID,"\t",$prot2,"\t",$range_list,"\t",$L,$extra,"\n";
}
close(SPEC);
exit(0);
