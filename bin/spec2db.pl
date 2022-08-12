#!/usr/bin/env perl

# Filename:         spec2db
# Description:      Exports DAIS ribosome spec/refs to a table format
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

use Digest::SHA qw(sha1_hex);
use strict;
use warnings;
use English qw( -no_match_vars );
use Carp qw(croak);

if ( scalar @ARGV != 2 ) {
    die("Usage:\n\tperl $PROGRAM_NAME <spec> <refs>\n\n");
}

sub nt_id2($) {
    my ($seq) = @_;

    if ( defined $seq ) {
        $seq =~ tr/ :.~-//d;
        if ( $seq ne q{} ) {
            return ( sha1_hex($seq), $seq );
        }
    } else {

        # Null string values for HIVE
        return ( '\N', '\N' );
    }
}

my %seqByRefSeg = ();
local $RS = '>';
my $FASTA;
open( $FASTA, '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");
while ( my $fasta_record = <$FASTA> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header   = shift(@lines);
    my $sequence = uc( join( q{}, @lines ) );

    my ( $refID, $seg, $id ) = ( split( '\|', $header ) )[0 .. 2];

    # exclude specially modified sequences
    if ( $id =~ /^\d+[A-Za-z]/smx || length($sequence) == 0 ) {
        next;
    }

    if ( $refID ne q{} && $seg ne q{} && !defined $seqByRefSeg{$refID}{$seg} ) {
        $seqByRefSeg{$refID}{$seg} = $sequence;
    }
}
close $FASTA or croak("Could not close $ARGV[1]: $OS_ERROR");

local $RS = "\n";
my $extra = q{};
my $SPEC;
open( $SPEC, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( my $line = <$SPEC> ) {
    chomp($line);
    my ( $ctype, $prot, $refID, $prot2, $range_list ) = split( "\t|[|]", $line );
    my $L      = 0;
    my @ranges = split( ';', $range_list );
    foreach my $pair (@ranges) {
        my ( $from, $to ) = split( ',', $pair );
        $L += $to - $from + 1;
    }
    $range_list =~ s/,/../gsmx;
    $range_list =~ tr/;/,/;

    if ( length( $seqByRefSeg{$refID}{$ctype} ) > 0 ) {
        $extra = "\t" . join( "\t", nt_id2( $seqByRefSeg{$refID}{$ctype} ) );
    } else {
        $extra = q{};
    }
    print STDOUT $ctype, "\t", $refID, "\t", $prot2, "\t", $range_list, "\t", $L, $extra, "\n";
}
close $SPEC or croak("Could not close $ARGV[0]: $OS_ERROR\n");
