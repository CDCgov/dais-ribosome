#!/usr/bin/env perl

# Filename:         makeInsProducts
# Description:      Creates CDS insertion table for aligned products from genome
#                   insertion table.
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

use English qw( -no_match_vars );
use strict;
use warnings;
use File::Basename qw(basename);
use File::Basename;
use Getopt::Long;
use Carp qw(croak);

my ( $geneSegment, $referenceID );
GetOptions( 'gene-segment|G=s' => \$geneSegment,
            'reference-id|R=s' => \$referenceID );

if ( scalar @ARGV != 3 ) {
    die("Usage:\n\t$PROGRAM_NAME <ins-table> <product-table> <prefix> [-G|--gene-segment]\n");
}

my $PROGRAM_SHORT_NAME = basename($PROGRAM_NAME);
my ( $insertionTable, $productTable, $prefix ) = @ARGV;
my ( $max, $productsFound ) = ( 0, 0 );
my $max_product_name = q{};
my %exons            = ();
my %fields           = ();
my %pMax             = ();

local $RS = "\n";
open( my $PROD, '<', $productTable ) or die("Cannot open $productTable for reading.\n");
while ( my $line = <$PROD> ) {
    chomp($line);
    my ( $segment, $peptide, $headerInfo, $coords ) = split( "\t", $line );

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
    $fields{$peptide} = $headerInfo;
    $pMax{$peptide}   = 0;

    foreach my $i ( 0 .. $#coordList ) {
        my ( $start, $stop ) = split( ',', $coordList[$i] );

        if ( $max < $stop ) {
            $max              = $stop;
            $max_product_name = $peptide;
        }
        if ( $pMax{$peptide} < $stop ) { $pMax{$peptide} = $stop; }
        $exons{$peptide}[$i][0] = $start;
        $exons{$peptide}[$i][1] = $stop;
    }
}
close $PROD or croak("Cannot close file $productTable: $OS_ERROR\n");

if ( defined $geneSegment && $productsFound == 0 ) {
    die("$PROGRAM_NAME:\tNo products found for $geneSegment.\n");
}

my @peptides = sort { $a cmp $b } keys(%fields);
my $filename = $prefix . '.ins';

local $RS = "\n";
open( my $OUT, '>', $filename )       or die("$PROGRAM_NAME:\tCannot open $filename for writing.\n");
open( my $INS, '<', $insertionTable ) or die("$PROGRAM_NAME:\tCannot open $insertionTable for reading.\n");

my @lines = <$INS>;
chomp(@lines);

foreach my $line (@lines) {
    my @fields = split( "\t", $line );
    if ( scalar @fields != 3 ) {
        die("$PROGRAM_NAME:\tExpected 3 fields in this format:\n\tID<tab>POSITION<tab>INSERT\n");
    } else {
        my ( $id, $pos, $insert ) = @fields;
        if ( $pos > $max ) {
            print STDERR "$PROGRAM_SHORT_NAME INFO:",
              "  Insertion exceeds range of annotated loci ($id): $pos >= $max ($max_product_name)\n";
            next;
        }

        foreach my $p (@peptides) {
            my $offset = 0;
            my $newPos = 0;
            foreach my $i ( 0 .. scalar( @{ $exons{$p} } ) - 1 ) {
                my ( $start, $stop ) = @{ $exons{$p}[$i] };
                if ( $start <= $pos && $pos < $stop ) {
                    $newPos = $pos - $start + 1 + $offset;
                    print $OUT $id, '|', $p, "\t", $newPos, "\t", $insert, "\n";
                    $offset += ( $stop - $start + 1 );
                    last;
                }
                $offset += ( $stop - $start + 1 );
            }
            if ( $pos == $max && $pos == $pMax{$p} ) {
                $newPos = $offset;
                print $OUT $id, '|', $p, "\t", $newPos, "\t", $insert, "\n";
            }
        }
    }
}
close $INS or croak("$PROGRAM_NAME:\tCannot close $insertionTable: $OS_ERROR\n");
close $OUT or croak("$PROGRAM_NAME:\tCannot close $filename: $OS_ERROR\n");
