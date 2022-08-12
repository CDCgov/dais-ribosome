#!/usr/bin/env perl

# Filename:         winnowSAM
# Description:      Pass thru the best scoring/matching query records with respect
#                   to multiple alignments to same-length references.
#
# Date dedicated:   2022-07-20
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Unpublished
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.

use strict;
use warnings;
use English qw( -no_match_vars );
use Getopt::Long;
use Carp qw(croak);

my ( $useMatches, $inPlace, $interleavedPairs );
GetOptions(
            'use-matches|M'       => \$useMatches,
            'in-place|I'          => \$inPlace,
            'interleaved-pairs|P' => \$interleavedPairs
);

if ( scalar @ARGV != 1 ) {
    die("\nUsage:\n\tperl $PROGRAM_NAME <sam>\n\n");
}

# FUNCTIONS #
sub countMatch($) {
    my ($cig) = @_;
    my $count = 0;
    while ( $cig =~ /(\d+)([MIDNSHP])/gsmx ) {
        if ( $2 eq 'M' ) {
            $count += $1;
        }
    }
    return $count;
}

#############
my $pair          = 0;
my %scoreByQuery  = ();
my %recordByQuery = ();
my $RS            = "\n";
my ( $previous, $header ) = ( q{}, q{} );

my $SAM;
open( $SAM, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( my $line = <$SAM> ) {
    if ( substr( $line, 0, 1 ) eq '@' ) {
        if ( $line ne $previous ) {
            $header .= $line;
        }
        $previous = $line;
        next;
    }

    chomp($line);
    my ( $qname, $flag, $rn, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual, $AS ) = split( "\t", $line );

    if ($interleavedPairs) {
        $qname = $qname . '_' . ( $pair % 2 );
        $pair++;
    }

    if ( $cigar eq '*' ) {
        next;
    }

    my $score;
    if ($useMatches) {
        $score = countMatch($cigar);
    } elsif ( $AS =~ /AS:\w:(\d+)/smx ) {
        $score = $1;
    } else {
        print STDERR "Warning, using simple matches as back-up: $qname\n";
        $score = countMatch($cigar);
    }

    if ( !defined $scoreByQuery{$qname} || $scoreByQuery{$qname} < $score ) {
        $scoreByQuery{$qname}  = $score;
        $recordByQuery{$qname} = $line;
    }
}
close $SAM or croak("Cannot close file $ARGV[0]: $OS_ERROR\n");

if ($inPlace) {
    my $SAM_OUT;
    open( $SAM_OUT, '>', $ARGV[0] ) or die("Cannot open $ARGV[0] for writing.\n");
    print $SAM_OUT $header;
    foreach my $query ( keys(%recordByQuery) ) {
        print $SAM_OUT $recordByQuery{$query}, "\n";
    }
    close $SAM_OUT or croak("Cannot close file $ARGV[0]: $OS_ERROR\n");
} else {
    print STDOUT $header;
    foreach my $query ( keys(%recordByQuery) ) {
        print STDOUT $recordByQuery{$query}, "\n";
    }
}
