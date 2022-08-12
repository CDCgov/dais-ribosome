#!/usr/bin/env perl

# Filename:         leftjoin
# Description:      Allows for left joins with multiple tables on a shared key.
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

use warnings;
use strict;
use English qw( -no_match_vars );
use Carp qw(croak);
use Getopt::Long;

my ( $delim, $fieldSet );
GetOptions( 'delim|D=s' => \$delim,
            'field|F=s' => \$fieldSet );

if ( scalar(@ARGV) < 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <table> [<left1> <left2> ...]\n"
         . "\t\t--delim|-D <CHAR>\tDelimiter for the key, the column delimiter many only be tab.\n"
         . "\t\t--field|-F <STR>\tComma-delimited set of fields to use for group. Default: column 1.\n\n" );
}

sub complementArray($$) {
    my ( $A1, $A2 ) = @_;
    my %H2  = map { $_ => 1 } @{$A2};
    my $key = q{};

    my @indices = grep { !defined( $H2{$_} ) } ( 0 .. $#{$A1} );
    return ( @{$A1}[@indices] );
}

my $maxSelected    = 1;
my $numberSelected = 0;
my @fields         = ();
if ( defined $fieldSet ) {
    @fields         = split( ',', $fieldSet );
    $numberSelected = scalar(@fields);
    foreach my $x (@fields) {
        if ( $x > $maxSelected ) { $maxSelected = $x; }
        if ( $x == 0 ) {
            die("$PROGRAM_NAME ERROR: field must be specified.\n");
        } elsif ( $x < 0 ) {
            die("$PROGRAM_NAME ERROR: field must be a positive number.\n");
        }
    }
    for my $x ( 0 .. $numberSelected - 1 ) { $fields[$x]--; }
} else {
    $fields[0] = 0;
}

if ( !defined $delim ) {
    $delim = '|';
} elsif ( $delim eq q{} ) {
    die("$PROGRAM_NAME ERROR: No delimiter argument detected.\n");
} elsif ( length($delim) > 1 ) {
    die("$PROGRAM_NAME ERROR: single character delimiter expected instead of '$delim'.\n");
}

my $numberFiles     = scalar(@ARGV);
my @data            = ();
my @lengthRemaining = ();
foreach my $i ( 1 .. $numberFiles - 1 ) {
    my $IN;
    open( $IN, '<', $ARGV[$i] ) or die("Cannot open file $ARGV[$i].\n");
    while ( my $line = <$IN> ) {
        chomp($line);
        my @values      = split( "\t", $line );
        my $numberFound = scalar(@values);
        my $id          = q{};
        if ( $numberSelected > 0 ) {
            if ( $maxSelected > $numberSelected ) {
                die(   "$PROGRAM_NAME ERROR: non-existant field specified."
                     . " Wanted $numberSelected (max: $maxSelected) but found $numberFound\n" );
            }
            $id = join( $delim, ( @values[@fields] ) );
        } else {
            $id = $values[$fields[0]];
        }

        if ( $id ne q{} ) {
            my @remainingColumns = map { $_ eq q{} ? '\N' : $_ } complementArray( \@values, \@fields );
            my $N                = scalar(@remainingColumns);
            if ( !defined $lengthRemaining[$i - 1] || $N > $lengthRemaining[$i - 1] ) {
                $lengthRemaining[$i - 1] = $N;
            }
            $data[$i - 1]{$id} = [@remainingColumns];
        }
    }
    close $IN or croak("Cannot close file $ARGV[$i]: $OS_ERROR\n");
}

my $IN;
open( $IN, '<', $ARGV[0] ) or die("Cannot open main table: $ARGV[0].\n");
while ( my $line = <$IN> ) {
    chomp($line);
    my @values      = split( "\t", $line );
    my $numberFound = scalar @values;
    my $id;
    if ( $numberSelected > 0 ) {
        if ( $maxSelected > $numberSelected ) {
            die(   "$PROGRAM_NAME ERROR: non-existant field specified."
                 . "Wanted $numberSelected (max: $maxSelected) but found $numberFound\n" );
        }
        $id = join( $delim, ( @values[@fields] ) );
    } else {
        $id = $values[$fields[0]];
    }

    if ( $id ne q{} ) {
        foreach my $i ( 1 .. $numberFiles - 1 ) {
            if ( defined $data[$i - 1]{$id} ) {
                my $N = scalar( @{ $data[$i - 1]{$id} } );
                if ( $N < $lengthRemaining[$i - 1] ) {
                    if ( $N > 0 ) {
                        $line .= "\t" . join( "\t", @{ $data[$i - 1]{$id} } );
                    }

                    $line .= "\t\\N" x ( $lengthRemaining[$i - 1] - $N );

                } else {
                    $line .= "\t" . join( "\t", @{ $data[$i - 1]{$id} } );
                }
            } else {
                $line .= "\t\\N" x ( $lengthRemaining[$i - 1] );
            }
        }
    }
    print STDOUT $line, "\n";
}
close $IN or croak("Cannot close $ARGV[0]: $OS_ERROR\n");
