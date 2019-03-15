#!/usr/bin/env perl
# SSS - 2019


$/ = "\n";
$L = <>;
chomp($L);

$ID = '\w+';
$annot = '[ABC](_[A-Z0-9]+){1,2}';
$seq = '[a-zA-Z.~-]+';

if ( $L =~ /^$ID\t$annot\t$seq$/ ) {
	print "atxt";
} elsif ( $L =~ /^$ID\t$seq$/ ) {
	print 'txt';
} elsif ( $L =~ /^>$ID\|$annot(\Z|\|)/ ) {
	print 'afa';
} elsif ( $L =~ /^>$ID$/ ) {
	print 'fa';
} else {
	print 'unk';
}

print "\n";
