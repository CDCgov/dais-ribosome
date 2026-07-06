#!/usr/bin/env perl
## no critic (RequireExtendedFormatting, RequireDotMatchAnything)
# In general:
# - Extended formatting has comments and what we have is readable enough.
# - Our use of `.` matching should not consume newlines.

use English qw(-no_match_vars);
use Carp    qw(croak);
use strict;
use warnings;

my $release_pkg = $ENV{RELEASE_PKG} // 'dais-ribosome-cli';

local $RS = undef;
open my $fh, '<', 'CHANGELOG.md' or croak "Can't open CHANGELOG.md: $OS_ERROR";
my $changelog = <$fh>;
close $fh or croak "Can't close CHANGELOG.md: $OS_ERROR";

my ( $version, $date );
if ( $changelog =~ /^## \[(\S*?)\] - (\S+?)$/m ) {
    ( $version, $date ) = ( $1, $2 );
} else {
    die "Could not find top changelog release heading like '## [x.y.z] - yyyy-mm-dd'\n";
}

local $RS = "\n";
my $pkgid = qx(cargo pkgid -p $release_pkg 2>&1);
if ( $CHILD_ERROR != 0 ) {
    die "cargo pkgid failed for package '$release_pkg': $pkgid\n";
}

my $toml_version;
if ( $pkgid =~ /[#@]([^@# \n]+?)$/m ) {
    $toml_version = $1;
} else {
    die "Could not parse version from cargo pkgid output: $pkgid\n";
}

if ( $date ne 'TBD' ) {
    if ( $date !~ /^20\d{2}-\d{2}-\d{2}$/m ) {
        die "Version $version has invalid date format: $date (expected yyyy-mm-dd or TBD)\n";
    }

    if ( $version ne $toml_version ) {
        die "Cargo.toml ($toml_version) mismatches changelog ($version / $date)!\n";
    }

    # Here we do need to consume newlines
    if ( $changelog !~ /<!-- Versions -->.*?^\[\Q$version\E\]:/sm ) {
        die "Version $version not linked!\n";
    }
} elsif ( $toml_version !~ /dev|rc\d+$/m ) {
    die "Cargo.toml version ($toml_version) should have a '-dev' suffix since the changelog is: ($version / $date)!\n";
}
