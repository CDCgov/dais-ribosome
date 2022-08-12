#!/usr/bin/env perl

# Filename:         amendInsertionTable
# Description:      Takes CDS insertion table and translates with appropriate
#                   codon numbers added to the table.
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
use Carp qw(croak);
use Getopt::Long;

my $filterPrefix;
my ( $useSTDOUT, $splitID ) = ( 0, 0 );
GetOptions(
            'use-std-out|U'         => \$useSTDOUT,
            'split-id-field|T'      => \$splitID,
            'apply-all-filters|F=s' => \$filterPrefix
);

if ( scalar @ARGV != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <insertion_table.txt>"
         . "[-A|--apply-all-filters <prefix>] [-U|--use-std-out] [-T|--split-id-field]\n\n" );
}

my $insertionTable = $ARGV[0];
my $RS             = "\n";
my %inserts        = ();

open( my $INS_IN, '<', $insertionTable ) or die("Cannot open $insertionTable for reading.\n");
my @lines = <$INS_IN>;
chomp(@lines);
foreach my $line (@lines) {
    my @fields = split( "\t", $line );
    if ( scalar @fields < 3 || scalar @fields > 4 ) {
        die("Expected 3 fields in this format:\n\tID<tab>POSITION<tab>INSERT\n");
    } elsif ( scalar @fields == 4 ) {
        die("Likely already translated.\n");
    } else {
        my ( $id, $pos, $insert ) = @fields;
        $inserts{$id}{$pos} = lc($insert);
    }
}
close $INS_IN or croak("Cannot close file $insertionTable: $OS_ERROR\n");

#<<< augmented translation table
my %gc = (
    'TAA'=>'*','TAG'=>'*','TAR'=>'*','TGA'=>'*','TRA'=>'*','GCA'=>'A','GCB'=>'A','GCC'=>'A','GCD'=>'A','GCG'=>'A','GCH'=>'A',
    'GCK'=>'A','GCM'=>'A','GCN'=>'A','GCR'=>'A','GCS'=>'A','GCT'=>'A','GCV'=>'A','GCW'=>'A','GCY'=>'A','TGC'=>'C','TGT'=>'C',
    'TGY'=>'C','GAC'=>'D','GAT'=>'D','GAY'=>'D','GAA'=>'E','GAG'=>'E','GAR'=>'E','TTC'=>'F','TTT'=>'F','TTY'=>'F','GGA'=>'G',
    'GGB'=>'G','GGC'=>'G','GGD'=>'G','GGG'=>'G','GGH'=>'G','GGK'=>'G','GGM'=>'G','GGN'=>'G','GGR'=>'G','GGS'=>'G','GGT'=>'G',
    'GGV'=>'G','GGW'=>'G','GGY'=>'G','CAC'=>'H','CAT'=>'H','CAY'=>'H','ATA'=>'I','ATC'=>'I','ATH'=>'I','ATM'=>'I','ATT'=>'I',
    'ATW'=>'I','ATY'=>'I','AAA'=>'K','AAG'=>'K','AAR'=>'K','CTA'=>'L','CTB'=>'L','CTC'=>'L','CTD'=>'L','CTG'=>'L','CTH'=>'L',
    'CTK'=>'L','CTM'=>'L','CTN'=>'L','CTR'=>'L','CTS'=>'L','CTT'=>'L','CTV'=>'L','CTW'=>'L','CTY'=>'L','TTA'=>'L','TTG'=>'L',
    'TTR'=>'L','YTA'=>'L','YTG'=>'L','YTR'=>'L','ATG'=>'M','AAC'=>'N','AAT'=>'N','AAY'=>'N','CCA'=>'P','CCB'=>'P','CCC'=>'P',
    'CCD'=>'P','CCG'=>'P','CCH'=>'P','CCK'=>'P','CCM'=>'P','CCN'=>'P','CCR'=>'P','CCS'=>'P','CCT'=>'P','CCV'=>'P','CCW'=>'P',
    'CCY'=>'P','CAA'=>'Q','CAG'=>'Q','CAR'=>'Q','AGA'=>'R','AGG'=>'R','AGR'=>'R','CGA'=>'R','CGB'=>'R','CGC'=>'R','CGD'=>'R',
    'CGG'=>'R','CGH'=>'R','CGK'=>'R','CGM'=>'R','CGN'=>'R','CGR'=>'R','CGS'=>'R','CGT'=>'R','CGV'=>'R','CGW'=>'R','CGY'=>'R',
    'MGA'=>'R','MGG'=>'R','MGR'=>'R','AGC'=>'S','AGT'=>'S','AGY'=>'S','TCA'=>'S','TCB'=>'S','TCC'=>'S','TCD'=>'S','TCG'=>'S',
    'TCH'=>'S','TCK'=>'S','TCM'=>'S','TCN'=>'S','TCR'=>'S','TCS'=>'S','TCT'=>'S','TCV'=>'S','TCW'=>'S','TCY'=>'S','ACA'=>'T',
    'ACB'=>'T','ACC'=>'T','ACD'=>'T','ACG'=>'T','ACH'=>'T','ACK'=>'T','ACM'=>'T','ACN'=>'T','ACR'=>'T','ACS'=>'T','ACT'=>'T',
    'ACV'=>'T','ACW'=>'T','ACY'=>'T','GTA'=>'V','GTB'=>'V','GTC'=>'V','GTD'=>'V','GTG'=>'V','GTH'=>'V','GTK'=>'V','GTM'=>'V',
    'GTN'=>'V','GTR'=>'V','GTS'=>'V','GTT'=>'V','GTV'=>'V','GTW'=>'V','GTY'=>'V','TGG'=>'W','TAC'=>'Y','TAT'=>'Y','TAY'=>'Y'
);
#>>>

my ( $REJ, $INS_OUT );
if ( defined $filterPrefix ) {
    open( $REJ,     '>', $filterPrefix . ".ins.filtered" ) or die("Cannot write ${filterPrefix}.ins.filtered\n");
    open( $INS_OUT, '>', $filterPrefix . ".ins" )          or die("Cannot write ${filterPrefix}.ins\n");
} elsif ($useSTDOUT) {
    $INS_IN = *STDOUT;
} else {
    open( $INS_OUT, '>', $insertionTable ) or die("Cannot open $insertionTable for writing.\n");
}

foreach my $id ( sort { $a cmp $b } keys(%inserts) ) {
    foreach my $nt_pos ( sort { $a <=> $b } keys( %{ $inserts{$id} } ) ) {
        my $insert = $inserts{$id}{$nt_pos};
        my $L      = length($insert);
        my $aa     = q{};
        if ( $L >= 3 ) {
            for ( my $i = 0; $i < $L; $i += 3 ) {    ## no critic (ControlStructures::ProhibitCStyleForLoops)
                my $codon = uc( substr( $insert, $i, 3 ) );
                if ( defined $gc{$codon} ) {
                    $aa .= $gc{$codon};
                } elsif ( length($codon) < 3 ) {
                    $aa .= '~';
                } else {
                    $aa .= 'X';
                }
            }
        } else {
            $aa = '?';
        }

        my $aa_pos   = int( $nt_pos / 3 );
        my $nt_shift = ( $nt_pos % 3 );

        if ($splitID) { $id = join( "\t", split( '\|', $id ) ); }

        if ( defined $filterPrefix ) {
            if ( $aa eq '?' || $insert =~ /^[nN]+$/smx ) {
                print $REJ $id, "\t", $aa_pos, "\t", $insert, "\t", $aa, "\t", $nt_pos, "\t", $nt_shift, "\n";
            } else {
                print $INS_OUT $id, "\t", $aa_pos, "\t", $insert, "\t", $aa, "\t", $nt_pos, "\t", $nt_shift, "\n";
            }
        } else {
            print $INS_OUT $id, "\t", $aa_pos, "\t", $insert, "\t", $aa, "\t", $nt_pos, "\t", $nt_shift, "\n";
        }
    }
}

if ( defined $filterPrefix ) {
    close $REJ or croak("Cannot close file: $OS_ERROR\n");
}
close $INS_OUT or croak("Cannot close file: $OS_ERROR\n");
