#!/usr/bin/env perl

# Filename:         makeProducts
# Description:      Splices nucleic acid CDS from aligned reference using product
#                   specification table.
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
use File::Basename;
use Getopt::Long;
use Carp qw(croak);

my ( $geneSegment, $referenceID );
GetOptions( 'gene-segment|G=s' => \$geneSegment,
            'reference-id|R=s' => \$referenceID );

if ( scalar @ARGV != 3 ) {
    die("Usage:\n\t$PROGRAM_NAME <fasta> <product-table> <prefix> [-G|--gene-segment <CTYPE>] [-R|--reference-id <STR>\n\n");
}

local $RS = "\n";
my $max           = 0;
my $productsFound = 0;
my %fields        = ();
my %context       = ();
my %exons         = ();

open( my $PROD, '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");
while ( my $line = <$PROD> ) {
    chomp($line);
    my ( $segment, $peptide, $headerInfo, $coords, $prefix, $suffix ) = split( "\t", $line );

    if ( !defined $prefix ) { $prefix = q{} }
    if ( !defined $suffix ) { $suffix = q{} }

    my $refID;
    if ( $headerInfo =~ /\|/smx ) {
        ($refID) = split( /\|/smx, $headerInfo );
    } else {
        $refID = $headerInfo;
    }

    if ( ( defined $geneSegment && $segment ne $geneSegment ) || ( defined $referenceID && $refID ne $referenceID ) ) {
        next;
    } else {
        $productsFound++;
    }

    my @coordList = split( ';', $coords );
    $fields{$peptide}  = $headerInfo;
    $context{$peptide} = [lc($prefix), lc($suffix)];
    foreach my $i ( 0 .. $#coordList ) {
        my ( $start, $stop ) = split( ',', $coordList[$i] );
        my $index = $start - 1;
        my $L     = $stop - $index;
        if ( $stop > $max ) {
            $max = $stop;
        }
        $exons{$peptide}[$i][0] = $index;
        $exons{$peptide}[$i][1] = $L;
    }
}
close $PROD or croak("Cannot close file $ARGV[1]: $OS_ERROR\n");

if ( defined $geneSegment && $productsFound == 0 ) {
    die("$PROGRAM_NAME:\tNo products found for $geneSegment.\n");
}

my @peptides = sort { $a cmp $b } keys(%fields);
my $filename = $ARGV[2] . '.products';

local $RS = '>';
open( my $OUT,   '>', $filename ) or die("Cannot open $filename for writing.\n");
open( my $FASTA, '<', $ARGV[0] )  or die("Cannot open $ARGV[0] for reading.\n");
while ( my $fasta_record = <$FASTA> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $id       = shift(@lines);
    my $sequence = lc( join( q{}, @lines ) );
    my $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    } elsif ( $max > $length ) {
        die("Found a sequence shorter ($length) than the last coordinate position ($max).\n");
    } else {
        foreach my $p (@peptides) {
            my ( $prefix, $suffix ) = @{ $context{$p} };
            my $cds = q{};
            foreach my $i ( 0 .. scalar( @{ $exons{$p} } ) - 1 ) {
                $cds .= substr( $sequence, $exons{$p}[$i][0], $exons{$p}[$i][1] );
            }

            if ( $prefix ne q{} && $prefix ne substr( $sequence, 0, length($prefix) ) ) {
                next;
            } elsif ( $suffix ne q{} && $suffix ne substr( $sequence, -length($suffix) ) ) {
                next;
            }

            $length = length($cds);
            if ( $length % 3 != 0 ) { die("$PROGRAM_NAME:\tNot in triplets ($length) for peptide '$p'.\n"); }
            print $OUT '>', $id, '|', $p, "\n", $cds, "\n";
        }
    }
}
close $FASTA or croak("Cannot close file: $OS_ERROR");
close $OUT   or croak("Cannot close file: $OS_ERROR");
