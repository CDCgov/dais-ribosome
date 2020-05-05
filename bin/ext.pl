#!/usr/bin/env perl
# SSS - 2019



if ( scalar(@ARGV) != 2 ) {
	die("$0 <input> <MODULE>\n\n");
}

open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].\n");
$/ = "\n";
my $module = $ARGV[1];
my $L = <IN>; chomp($L);
close(IN);

my $ID = '\w+';
my $annot = '[ABC](_[A-Z0-9]+){1,2}';
my $seq = '[a-zA-Z.~-]+';

if ( $module =~ /CORONAVIRUS/i ) {
	$annot = '[A-Z]+-CoV(-\w+)*';
	$ID = '[A-Za-z0-9_-]+';
}

if ( $L =~ /^$ID\t$annot\t$seq$/ ) {
	print "atxt";
} elsif ( $L =~ /^$ID\t$seq$/ ) {
	print 'txt';
} elsif ( $L =~ /^>$ID\|$annot(\r?\z|\|)/ ) {
	print 'afa';
} elsif ( $L =~ /^>$ID\r?$/ ) {
	print 'fa';
} else {
	print 'unk';
}

print "\n";
