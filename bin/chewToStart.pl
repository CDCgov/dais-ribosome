#!/usr/bin/env perl

# Filename:         chewToStart
# Description:      For query sequences longer than reference with start codon
#                   trim the query sequence to the start codon.
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
#  material

use strict;
use warnings;
use English qw( -no_match_vars );

my ( $referenceFile, $referenceLength, $flexLen );
use Getopt::Long;
GetOptions(
            'ref|R=s'         => \$referenceFile,
            'length|L=i'      => \$referenceLength,
            'flex-length|X=i' => \$flexLen
);

if ( -t STDIN && scalar(@ARGV) != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <nts.fasta> [options]\n"
         . "\t\t-R|--ref\t\tReference file name.\n"
         . "\t\t-L|--length\t\tReference length.\n"
         . "\t\t-X|--flex-length\tAmount the 'chewed' sequence may be less than the reference.\n\n" );
}

my $REF_LEN = 0;
local $RS = '>';
if ( defined $referenceFile ) {
    my $REF;
    open( $REF, '<', $referenceFile ) or die("Cannot open $referenceFile for reading.\n");
    while ( my $fasta_record = <$REF> ) {
        chomp($fasta_record);
        my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
        my $id       = shift(@lines);
        my $sequence = uc( join( q{}, @lines ) );
        my $length   = length($sequence);
        if ( $length == 0 ) {
            next;
        } else {
            $REF_LEN = $length;
            last;
        }
    }
} elsif ( defined $referenceLength && $referenceLength > 0 ) {
    $REF_LEN = $referenceLength;
} else {
    $REF_LEN = 0;
}

if ( defined $flexLen ) { $REF_LEN -= $flexLen; }
while ( my $fasta_record = <> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $id       = shift(@lines);
    my $sequence = uc( join( q{}, @lines ) );
    my $length   = length($sequence);
    if ( $length == 0 ) { next; }

    if ( $length > $REF_LEN && $sequence =~ /ATG/smx ) {
        my $newSeq = substr( $sequence, $LAST_MATCH_START[0] );
        my $newLen = length($newSeq);
        if ( $newLen >= $REF_LEN ) {
            print STDOUT '>', $id, "\n", $newSeq, "\n";
        } else {
            print STDOUT '>', $id, "\n", $sequence, "\n";
        }
    } else {
        print STDOUT '>', $id, "\n", $sequence, "\n";
    }
}
