#!/usr/bin/env perl
use English qw(-no_match_vars);
use Carp    qw(croak);
use strict;
use warnings;

my $release_pkg = $ENV{RELEASE_PKG} // 'dais-ribosome-cli';

local $RS = undef;
open my $fh, '<', 'CHANGELOG.md' or croak "Can't open CHANGELOG.md: $OS_ERROR";
my $changelog = <$fh>;
close $fh or croak "Can't close CHANGELOG.md: $OS_ERROR";

my ( $version, $date ) = $changelog =~ /^## \[(.*?)\] - (\S+?)$/m
  or die "Could not find top changelog release heading like '## [x.y.z] - yyyy-mm-dd'\n";

local $RS = "\n";
my $pkgid = qx(cargo pkgid -p $release_pkg 2>&1);
die "cargo pkgid failed for package '$release_pkg': $pkgid" if $CHILD_ERROR != 0;

my ($toml_version) = $pkgid =~ /[#@]([^ \n]+)$/;
die "Could not parse version from cargo pkgid output: $pkgid" if !defined $toml_version;

if ( $date ne 'TBD' ) {
    if ( $date !~ /^20\d{2}-\d{2}-\d{2}$/smx ) {
        die "Version $version has invalid date format: $date (expected yyyy-mm-dd or TBD)\n";
    }

    if ( $version ne $toml_version ) {
        die "Cargo.toml ($toml_version) mismatches changelog ($version / $date)!\n";
    }
} elsif ( $toml_version !~ /dev|rc\d+$/smx ) {
    die "Cargo.toml version ($toml_version) should have a '-dev' suffix since the changelog is: ($version / $date)!\n";
}

if ( $changelog !~ /<!-- Versions -->.*?^\[\Q$version\E\]:/sm ) {
    die "Version $version not linked!\n";
}
