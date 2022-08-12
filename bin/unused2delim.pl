#!/usr/bin/env perl

# Filename:         unused2delim
# Description:      Converts fasta to a delimited format for unused data and null
#                   pads where appropriate.
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
use English qw( -no_match_vars);
use Getopt::Long;

#  NULL constant useful as special value
## no critic (ValuesAndExpressions::ProhibitConstantPragma)
use constant NULL => '\N';

my $fieldsExpected;
GetOptions( 'fields-expected|F=i' => \$fieldsExpected );

if ( -t STDIN && !scalar @ARGV ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <annotated.fasta> [options]\n"
         . "\t\t-F|--fields-expected <+INT>\t\tPads with nulls up to specified number of fields.\n\n" );
}

# Trim function.
# Removes whitespace from the start and end of the string
# TO-DO Use of uninitialized value $string in pattern match (m//) at /home/vfn4/dev/dais-ribosome/bin/unused2delim.pl line 59, <> chunk 1.
sub trim($) {
    my $string = shift;
    if ( defined $string && $string =~ /^\s*(.*?)\s*$/smx ) {
        return $1;
    } else {
        return $string;
    }
}

my @nullpad = ();
if ( defined $fieldsExpected && int $fieldsExpected > 0 ) {
    my $N = int $fieldsExpected;
    foreach my $i ( 0 .. $N ) {
        $nullpad[$i] = NULL;
    }
} else {
    $fieldsExpected = 1;
}
my $limit = $fieldsExpected - 1;

local $RS = '>';
while ( my $fasta_record = <> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header   = trim( shift(@lines) );
    my $sequence = uc( join( q{}, @lines ) );

    if ( length($sequence) == 0 ) { next; }

    my @fields = split( '\|', $header );
    my $N      = scalar(@fields);
    if ( $N < $fieldsExpected ) {
        my $diff = $fieldsExpected - $N - 1;
        print STDOUT join( "\t", @fields ), "\t", join( "\t", @nullpad[0 .. $diff] ), "\n";
    } else {
        print STDOUT join( "\t", @fields[0 .. $limit] ), "\n";
    }
}
