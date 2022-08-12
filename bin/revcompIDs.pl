#!/usr/bin/env perl

# Filename:         revcompID
# Description:      Applies reverse complement to [ID<TAB>...Strand] text file.
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

my ( $inPlace, $strandField );
GetOptions( 'in-place|I' => \$inPlace, 'strand-field|S' => \$strandField );

if ( scalar @ARGV != 2 ) {
    die(  "Usage:\n\tperl $PROGRAM_NAME <fasta> <IDs> [options]\n"
        . "\t\t-I|--in-place\t\tWrite out file in place rather than to STDOUT.\n"
        . "\t\t-S|--strand-field\tExpects tab-delimited fields (ID:first, strand:last). If strand is '-', reverse complement.\n\n"
    );
}

local $RS = "\n";
my %IDs = ();
my $TAB;
open( $TAB, '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");
while ( my $line = <$TAB> ) {
    chomp($line);
    my @f = split( "\t", $line );
    if ( defined $f[0] && $f[0] ne q{} ) {
        if ( !defined $strandField || $f[-1] eq '-' ) {
            $IDs{ $f[0] } = 1;
        }
    }
}
close $TAB or croak("Cannot close file $ARGV[1]: $OS_ERROR");

local $RS = '>';
if ( defined $inPlace ) {
    open( my $FASTA, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
    my @records = <$FASTA>;
    chomp(@records);
    close($FASTA) or croak("Cannot close file $ARGV[0]: $OS_ERROR\n");

    open( $FASTA, '>', $ARGV[0] ) or die("Cannot open $ARGV[0] for writing.\n");
    foreach my $fasta_record (@records) {
        my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
        my $seqID    = shift(@lines);
        my $sequence = lc( join( q{}, @lines ) );

        if ( length($sequence) == 0 ) { next; }

        if ( defined $IDs{$seqID} ) {
            $sequence = reverse($sequence);
            $sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
        }

        print $FASTA '>', $seqID, "\n", $sequence, "\n";
    }
} else {
    open( my $FASTA, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
    while ( my $fasta_record = <$FASTA> ) {
        chomp($fasta_record);
        my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
        my $seqID    = shift(@lines);
        my $sequence = lc( join( q{}, @lines ) );

        if ( length($sequence) == 0 ) {
            next;
        }

        if ( defined $IDs{$seqID} ) {
            $sequence = reverse($sequence);
            $sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
        }

        print STDOUT '>', $seqID, "\n", $sequence, "\n";
    }
    close $FASTA or croak("Cannot close file $ARGV[0]: $OS_ERROR\n");
}
