#!/usr/bin/env perl
# Samuel Shepard - 2019.03
# Annotate the CDS of sequences in DAIS

use warnings;
use strict;
use English qw(-no_match_vars);
use File::Basename;
use Carp qw(croak);

if ( scalar @ARGV != 6 ) {
    die(   "\nUsage:\n\t$PROGRAM_NAME <spec> <original-fasta>"
         . " <product-fasta-c> <product-ins-c>"
         . " <seg-alignment-fasta> <seg-insertion-ins-txt>"
         . "\n\n" );
}

my ( $specfile, $origfile, $fasprodfile, $insprodfile, $fassegfile, $inssegfile ) = (@ARGV);

# FUNCTIONS #
sub removeElongation($) {
    my ($seq) = @_;
    $seq =~ s/[A-Z]+$//smx;
    return $seq;
}

sub sequenceByCigar($$$) {
    my ( $sequence, $cigar, $offset ) = @_;
    my ( $result, $state ) = ( q{}, q{} );
    my $length = 0;
    while ( $cigar =~ m/(\d+)([MIDSHN])/gsmx ) {
        ( $length, $state ) = ( $1, $2 );
        if ( $state =~ /[MI]/smx ) {
            $result .= substr( $sequence, $offset, $length );
            $offset += $length;
        }
    }
    return $result;
}

sub getSubstringOffset($$) {    ## no critic (Subroutines::RequireArgUnpacking)
    my ( $original, $alignable ) = ( lc( $_[0] ), lc( $_[1] ) );
    my $leftpad = 0;
    if ( $alignable =~ /^(\.+)/smx ) {
        $leftpad = length($1);
    }
    $alignable =~ tr/.-//d;

    if ( $original =~ m/\Q$alignable\E/smx ) {
        return ( $LAST_MATCH_START[0] - $leftpad );
    } else {
        print STDERR "ERROR, alignable not found.\n> Original(gSO)\n$original\n\n>Alignable\n$alignable\n\n";
        return 0.5;
    }
}

sub getSubstringCoords($$) {    ## no critic (Subroutines::RequireArgUnpacking)
    my ( $original, $alignable ) = ( lc( $_[0] ), lc( $_[1] ) );
    $alignable =~ tr/.-//d;

    if ( $original =~ m/\Q$alignable\E/smx ) {
        return ( $LAST_MATCH_START[0], $LAST_MATCH_END[0] );
    } else {
        print STDERR "ERROR, alignable not found.\n>Original(gSC)\n$original\n\n>Alignable\n$alignable\n\n";
        return ( -1, -1 );
    }
}

sub sequenceToCigar($) {
    my ($seq) = @_;
    $seq =~ tr/A-Z/I/;
    $seq =~ tr/a-z/M/;
    $seq =~ tr/-/D/;
    $seq =~ tr/./N/;

    my $cigar = q{};
    while ( $seq =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/gsmx ) {
        $cigar .= length($1) . substr( $1, 0, 1 );
    }
    return $cigar;
}

sub condenseCigar($) {
    my ($cig) = @_;
    my $cigar = q{};
    my $state = q{};
    while ( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/gsmx ) {
        $state = $1;
        $cigar .= length($state);
        $cigar .= substr( $state, 0, 1 );
    }
    return $cigar;
}

sub sequenceToStates($) {
    my ($seq) = @_;
    $seq =~ tr/A-Z/I/;
    $seq =~ tr/a-z/M/;
    $seq =~ tr/-/D/;
    $seq =~ tr/./N/;
    return $seq;
}

sub addInsertions($$) {    ## no critic (Subroutines::RequireArgUnpacking)
    my $seq     = lc( $_[0] );
    my $inserts = $_[1];
    my ( $offset, $pos ) = 0;
    my $insert = q{};
    foreach my $pos ( sort { $a <=> $b } keys( %{$inserts} ) ) {
        $insert = $inserts->{$pos};
        substr( $seq, int($pos) + $offset, 0, uc($insert) );
        $offset += length($insert);
    }

    return $seq;
}

sub addInsertionsBounded($$$) {    ## no critic (Subroutines::RequireArgUnpacking)
    my $seq     = lc( $_[0] );
    my $inserts = $_[1];
    my $offset  = $_[2];
    my $pos     = 0;
    my $L       = length($seq);
    foreach my $pos ( sort { $a <=> $b } keys( %{$inserts} ) ) {

        # 1 - based check
        if ( ( $pos + $offset ) > $L || ( $pos + $offset ) < 1 ) { last; }
        substr( $seq, $pos + $offset, 0, uc( $inserts->{$pos} ) );
        $offset += length( $inserts->{$pos} );
        $L = length($seq);
    }

    return $seq;
}

# process specifications
local $RS = "\n";
my %specs           = ();
my %segmentByRefPep = ();
my ( $max, $productsFound ) = ( 0, 0 );

open( my $SPEC, '<', $specfile ) or croak("$PROGRAM_NAME ERROR: Could not open $specfile for reading.\n");
while ( my $line = <$SPEC> ) {
    chomp($line);
    my ( $segment, $peptide, $headerInfo, $coords, $prefix, $suffix ) = split( "\t", $line );
    my ( $ref_id, $peptide2 ) = split( '\|', $headerInfo );
    my @coordList = split( ';', $coords );

    $segmentByRefPep{$ref_id}{$peptide} = $segment;
    foreach my $i ( 0 .. $#coordList ) {
        my ( $start, $stop ) = split( ',', $coordList[$i] );
        $specs{$ref_id}{$peptide}[$i][0] = $start - 1;
        $specs{$ref_id}{$peptide}[$i][1] = $stop - $start + 1;
    }
}
close $SPEC or croak("Cannot close file: $OS_ERROR\n");

# process original fasta
local $RS = q{>};
my %originals = ();

open( my $ORIG, '<', $origfile ) or die("$PROGRAM_NAME ERROR: Could not open $origfile for reading.\n");
while ( my $fasta_record = <$ORIG> ) {
    chomp($fasta_record);
    my @lines  = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header = shift(@lines);
    my $seq    = lc( join( q{}, @lines ) );

    if ( length($seq) == 0 ) {
        next;
    } else {
        my @id = split( '\|', $header );
        $originals{ $id[0] } = $seq;
    }
}
close $ORIG or croak("Cannot close file $origfile: $OS_ERROR\n");

# process alignment (gene segment level) insertions
local $RS = "\n";
my %segmentInsertions = ();

open( my $SEG_INS, '<', $inssegfile ) or die("$PROGRAM_NAME ERROR: could not open $inssegfile for reading.\n");
while ( my $line = <$SEG_INS> ) {
    chomp($line);
    my ( $compound_id, $upstream_position, $insert ) = split( "\t", $line );
    my ( $flu_seq_id,  $segment,           $ref_id ) = split( '\|', $compound_id );
    $segmentInsertions{$ref_id}{$segment}{$flu_seq_id}{$upstream_position} = $insert;
}
close $SEG_INS or croak("Cannot close file $inssegfile: $OS_ERROR\n");

# [INPUT]: "240934|A_HA_H3|HK4801|PB1	2274	TGA"
# process product (nucleotide gene transcript) insertions
local $RS = "\n";
my %productInsertions = ();
open( my $PROD_INS, '<', $insprodfile ) or die("$PROGRAM_NAME ERROR: could not open $insprodfile for reading.\n");
while ( my $line = <$PROD_INS> ) {
    chomp($line);
    my ( $compound_id, $upstream_position, $insert ) = split( "\t", $line );
    my ( $flu_seq_id, $segment, $ref_id, $peptide ) = split( '\|', $compound_id );
    $productInsertions{$ref_id}{$peptide}{$flu_seq_id}{$upstream_position} = $insert;
}
close $PROD_INS or croak("Cannot close file $insprodfile: $OS_ERROR\n");

# process the segment coordinates to create a bounding box for segment alignments via reference
local $RS = ">";
my %segmentOffset = ();
my ( $total, $found ) = ( 0, 0 );
open( my $ALIGNED, '<', $fassegfile ) or die("$PROGRAM_NAME ERROR: Could not open $fassegfile for reading.\n");
while ( my $fasta_record = <$ALIGNED> ) {
    chomp($fasta_record);
    my @lines  = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header = shift(@lines);
    my $seq    = lc( join( q{}, @lines ) );

    if ( length($seq) == 0 ) { next; }

    my ( $flu_seq_id, $segment, $ref_id ) = split( '\|', $header );

    if ( defined $originals{$flu_seq_id} ) {
        $found++;

        # Add back in insertions so the reference-aligned sequence will map to the original.
        # Remove the trailing insertions (3' elongation) if applicable.
        # Goal is to get the Reference boundary coords, which excludes elongation.
        if ( defined $segmentInsertions{$ref_id}{$segment}{$flu_seq_id} ) {
            $seq = removeElongation( addInsertions( $seq, \%{ $segmentInsertions{$ref_id}{$segment}{$flu_seq_id} } ) );
        }

        # Reference-Aligned Query to Original
        my $offset = getSubstringOffset( $originals{$flu_seq_id}, $seq );
        if ( $offset == 0.5 ) {
            die("Issue with $segment / $ref_id / #$flu_seq_id!\n");
        }
        $segmentOffset{$segment}{$ref_id}{$flu_seq_id} = $offset;
    } else {
        print STDERR "Original pair not found: $flu_seq_id ( $segment / $ref_id )\n";
    }

    $total++;
}
close $ALIGNED or croak("Cannot close file fassegfile: $OS_ERROR\n");

# Sample header: >251324|A_NA_N2|HK4801|NA
# process final products and create a map between reference coordinates (codon number) and original nucleotide coordinates
local $RS = ">";
open( my $PRODUCTS, '<', $fasprodfile ) or die("$PROGRAM_NAME ERROR: Could not open $fasprodfile for reading.\n");
while ( my $fasta_record = <$PRODUCTS> ) {
    chomp($fasta_record);
    my @lines       = split( /\r\n|\n|\r/smx, $fasta_record );
    my $compound_id = shift(@lines);
    my $seq_prod    = lc( join( q{}, @lines ) );
    my $seq_len     = length($seq_prod);

    if ( $seq_len == 0 ) { next; }

    my ( $flu_seq_id, $segment, $ref_id, $peptide ) = split( '\|', $compound_id );
    my $leftPad = $seq_prod =~ /^(\.+)/smx ? length($1) : 0;

    if ( defined $originals{$flu_seq_id} ) {
        my @exons    = @{ $specs{$ref_id}{$peptide} };
        my $original = $originals{$flu_seq_id};

        # Aligned segment to original sequence offset. (1-based).
        my $oriOffset;
        if ( defined $segmentOffset{$segment}{$ref_id}{$flu_seq_id} ) {
            $oriOffset = $segmentOffset{$segment}{$ref_id}{$flu_seq_id};
        } else {
            die("Missing segment offset data for $segment / $ref_id / $flu_seq_id\n");
        }
        my ( $oriCoords, $pepCoords ) = ( q{}, q{} );
        my $pepOffset = 0;
        my $first     = 1;
        foreach my $i ( 0 .. $#exons ) {
            my ( $idx, $L ) = @{ $exons[$i] };

            # Get the cigar for the current exon, do not add insertions outside the valid bounds.
            # We process the exons from the peptide / product sequence.
            # I use the peptide term interchangably with product, although the product files are untranslated in DAIS terms.
            my $exon = substr( $seq_prod, $pepOffset, $L );
            my $exonCigar = sequenceToCigar(
                      addInsertionsBounded( $exon, \%{ $productInsertions{$ref_id}{$peptide}{$flu_seq_id} }, -$pepOffset ) );

            # Let $idx + 1 = offset for the peptide within the segment alignment
            # Thus the starting original coordinate adds the peptide to ref and ref to original offsets.
            my $oriCursor = $idx + 1 + $oriOffset + $leftPad;

            # Peptide cursor is after previous exon lengths (offsets). Start at 1.
            my $pepCursor = $pepOffset + 1;

            # Create a map from the original coordinate space to the peptide product coordinate space.
            # Based on the OPeration, we advance each coordinate system relative to the other.
            while ( $exonCigar =~ m/(\d+)([MDNI])/gsmx ) {
                my ( $inc, $op ) = ( $1, $2 );
                if ( $op eq 'N' ) {    ## no critic (ControlStructures::ProhibitCascadingIfElse)
                    if ($first) {
                        $pepCursor += $inc;
                        $first = 0;
                    }
                } elsif ( $op eq 'M' ) {
                    $pepCoords .= $pepCursor . '..';
                    $oriCoords .= $oriCursor . '..';

                    $pepCursor += $inc;
                    $oriCursor += $inc;

                    $pepCoords .= ( $pepCursor - 1 ) . ';';
                    $oriCoords .= ( $oriCursor - 1 ) . ';';
                } elsif ( $op eq 'I' ) {
                    $pepCoords .= ( $pepCursor - 1 ) . ';';
                    $oriCoords .= $oriCursor . '..';

                    $oriOffset += $inc;    # Take into account insertions in the original sequence.
                    $oriCursor += $inc;

                    $oriCoords .= ( $oriCursor - 1 ) . ';';
                } elsif ( $op eq 'D' ) {
                    $pepCursor += $inc;
                    $oriOffset -= $inc;    # Likewise deletions are removed from the offset.
                } else {
                    die("$op : Unknown\n");
                }
            }
            $pepOffset += $L;    # reflect that the exon has moved forward
        }

        chop($pepCoords);
        chop($oriCoords);
        print STDOUT join( '|', ( $flu_seq_id, $segment, $ref_id, $peptide ) ), "\t", $oriCoords, "\t", $pepCoords, "\n";
    } else {
        print STDERR "Original pair not found: $flu_seq_id ( $segment / $peptide / $ref_id )\n";
    }
}
close $PRODUCTS or croak("Cannot close file $fasprodfile: $OS_ERROR\n");
