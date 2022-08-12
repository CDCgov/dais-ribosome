#!/usr/bin/env perl

# Filename:         ext
# Description:      File input detection by printing extension code. The "a" is
#                   for annotated.
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
use English qw(-no_match_vars);
use Carp qw(croak);

if ( scalar @ARGV != 2 ) {
    die("$PROGRAM_NAME <input> <MODULE>\n\n");
}

my $IN;
open( $IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
local $RS = "\n";
my $module = $ARGV[1];
my $L      = <$IN>;
chomp($L);
close $IN or croak("Cannot close file: $OS_ERROR\n");

my $ID    = '\w+';
my $annot = '[ABC](_[A-Z0-9]+){1,2}';
my $seq   = '[a-zA-Z.~-]+';

if ( $module =~ /CORONAVIRUS/ismx ) {
    $annot = '[A-Z]+-CoV(-\w+)*';
    $ID    = '[A-Za-z0-9_-]+';
}

my $type = 'unk';
if ( $L =~ /^$ID\t$annot\t$seq$/smx ) {    ## no critic (ControlStructures::ProhibitCascadingIfElse)
    $type = "atxt";
} elsif ( $L =~ /^$ID\t$seq$/smx ) {
    $type = 'txt';
} elsif ( $L =~ /^>$ID\|$annot(\r?\Z|\|)/smx ) {
    $type = 'afa';
} elsif ( $L =~ /^>$ID\r?\Z/smx ) {
    $type = 'fa';
}

print STDOUT "$type\n";