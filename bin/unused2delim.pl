#!/usr/bin/env perl
# Converts fasta to a delimited format for unused data


use constant NULL => '\N';
use Getopt::Long;
GetOptions( 
		'fields-expected|F=i' => \$fieldsExpected
	);

if ( -t STDIN && ! scalar(@ARGV) ) {
	$message = "Usage:\n\tperl $0 <annotated.fasta> [options]\n";
	$message .= "\t\t-F|--fields-expected <+INT>\t\tPads with nulls up to specified number of fields.\n";
	die($message."\n");
}

my @nullpad = ();
if ( defined($fieldsExpected) && int($fieldsExpected) > 0 ) {
	$N = int($fieldsExpected);
	foreach my $i ( 0 .. $N ) {
		$nullpad[$i] = NULL;
	}
} else {
	$fieldsExpected = 1;
}
$limit = $fieldsExpected - 1;

$/ = '>';
while( $record = <> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = trim(shift(@lines));
	$sequence = uc(join('',@lines));

	if ( length($sequence) == 0 ) { next; }

	@fields = split('\|',$header);
	$N = scalar(@fields);
	if ( $N < $fieldsExpected ) {
		$diff = $fieldsExpected - $N -1;
		print STDOUT join("\t",@fields),"\t",join("\t",@nullpad[0..$diff]),"\n";
	} else {
		print STDOUT join("\t",@fields[0..$limit]),"\n";
	}
}

# Trim function.
# # Removes whitespace from the start and end of the string
sub trim($) {
 	my $string = shift;
	$string =~ /^\s*(.*?)\s*$/;
 	return $1;
}
