#!/usr/bin/env perl
# Filename:         filterUntranslatedInserts
# Description:      Filters inserts outside the translated range.
#
# Date dedicated:   2023-06-26
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
use English qw(-no_match_vars);
use Carp qw(croak);

if ( scalar( @ARGV != 2 ) ) {
    die("\n\tUsage: perl $PROGRAM_NAME <coords_file> <inserts_file>\n\n\tFilters on STDIN to STDOUT!\n");
}

local $RS = "\n";
my %last_position_by_id = ();
open( my $COORD, '<', $ARGV[0] ) or die("Cannot open for reading: $OS_ERROR\n");

while ( my $coord_line = <$COORD> ) {
    chomp($coord_line);
    my ( $id, $query_coords, $cds_coords ) = split( "\t", $coord_line );
    if ( $cds_coords =~ /\D(\d+)$/smx ) {
        $last_position_by_id{$id} = $1;
    }

}
close($COORD) or croak("Cannot close file: $OS_ERROR\n");

open( my $INSERTS, '<', $ARGV[1] ) or croak("Cannot open for reading: $OS_ERROR\n");
while ( my $insert_line = <$INSERTS> ) {
    chomp($insert_line);
    my ( $id, $residue, $insertedNT, $insertedAA, $query_pos, $frame ) = split( "\t", $insert_line );

    # We allow hanging inserts.
    if ( !defined $last_position_by_id{$id} || $query_pos <= $last_position_by_id{$id} ) {
        print $insert_line, "\n";
    }

}
close($INSERTS) or croak("Cannot close file: $OS_ERROR\n");