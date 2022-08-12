#!/usr/bin/env perl

# Filename:         seq2del
# Description:      Creates deletion table from aligned sequence file.
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

use POSIX;
use strict;
use warnings;
use English qw( -no_match_vars );
use Getopt::Long;
use Carp qw(croak);

my $genomeDeletions = 0;
GetOptions( 'genome|G' => \$genomeDeletions );

my $AA_ALN  = 6;
my $CDS_ID  = 7;
my $CDS_ALN = 11;
my $GEN_ALN = 7;

my $aln     = $CDS_ALN;
my $pre_end = 4;
if ($genomeDeletions) {
    $aln     = $GEN_ALN;
    $pre_end = 2;
}

if ( scalar @ARGV != 1 ) {
    die("\nUsage:\n\t$PROGRAM_NAME <input.seq> [--genome|-G]\n\n");
}

local $RS = "\n";

my $IN;
open( $IN, '<', $ARGV[0] ) or die("$PROGRAM_NAME ERROR: cannot open $ARGV[0] for reading!\n");
while ( my $line = <$IN> ) {
    chomp($line);
    my @f = split( "\t", $line );

    if ( !defined $f[$aln] ) {
        next;
    }

    my $prefix = join( "\t", @f[0 .. $pre_end] );
    while ( $f[$aln] =~ m/(-+)/gsmx ) {
        my ( $b, $e, $l ) = ( $LAST_MATCH_START[0] + 1, $LAST_MATCH_END[0], $LAST_MATCH_END[0] - $LAST_MATCH_START[0] );
        my $in_frame = $b % 3 == 1 && $e % 3 == 0 ? 'true' : 'false';
        my ( $aa_b, $aa_e, $aa_l ) = ( int( ( $b - 1 ) / 3 + 1 ), int( ( $e - 1 ) / 3 + 1 ), ceil( $l / 3 ) );

        if ($genomeDeletions) {
            print STDOUT $prefix, "\t", $b, "\t", $e, "\t", $l, "\n";
        } else {
            print STDOUT $prefix, "\t", $aa_b, "\t", $aa_e, "\t", $aa_l, "\t", $in_frame, "\t", $f[$CDS_ID], "\t", $b, "\t",
              $e, "\t", $l, "\n";
        }
    }
}
close($IN) or croak("File could not close: $OS_ERROR")
