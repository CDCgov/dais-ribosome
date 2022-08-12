#!/usr/bin/env perl

# Filename:         getSSWconfig
# Description:      Gets argument list for SSW for shell execution based on gene.
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

if ( scalar(@ARGV) != 2 ) {
    die("\nUsage:\n\tperl $PROGRAM_NAME <config_file> <gene_segment_subtype>\n\n");
}

my $filename      = $ARGV[0];
my $gene_selected = $ARGV[1];

local $RS = "\n";
open( my $IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
my $default = <$IN>;
chomp($default);
while ( my $line = <$IN> ) {
    chomp($line);
    my ( $gene, $params, $ref ) = split( "\t", $line );
    if ( index( $gene_selected, $gene ) > -1 ) {
        my ( $match, $mismatch, $gapopen, $gapextend ) = split( q{ }, $params );
        if ( $match < 1 )     { $match     = 1; }
        if ( $mismatch < 1 )  { $mismatch  = 1; }
        if ( $gapopen < 1 )   { $gapopen   = 1; }
        if ( $gapextend < 1 ) { $gapextend = 1; }

        print " -m $match -x $mismatch -o $gapopen -e $gapextend ";
        exit;
    }
}

my ( $gene, $params, $ref ) = split( "\t", $default );
my ( $match, $mismatch, $gapopen, $gapextend ) = split( q{ }, $params );
if ( $match < 1 )     { $match     = 1; }
if ( $mismatch < 1 )  { $mismatch  = 1; }
if ( $gapopen < 1 )   { $gapopen   = 1; }
if ( $gapextend < 1 ) { $gapextend = 1; }

print STDOUT " -m $match -x $mismatch -o $gapopen -e $gapextend ";

