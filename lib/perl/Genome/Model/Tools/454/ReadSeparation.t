#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;


BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    }
    plan skip_all => "skipping this test for now";
    use_ok('Genome::Model::Tools::454::ReadSeparation');
    use_ok('Genome::Model::Tools::454::IsolatePrimerTag');
    use_ok('Genome::Model::Tools::454::CrossMatchPrimerTag');
    use_ok('Genome::Model::Tools::454::SeparateReadsWithCrossMatchAlignment');
}

my $sff_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-454-ReadSeparation/test_454_primer_tag_100k.sff';

my @version_subdirs = qw/ offInstrumentApps mapasm454_source /;
foreach my $sub_dir (@version_subdirs) {
    my $version;
    $version = '2.0.00.20-1' if $sub_dir eq 'offInstrumentApps';
    $version = '10282008' if $sub_dir eq 'mapasm454_source';
    my $read_separation = Genome::Model::Tools::454::ReadSeparation->create(
									    sff_file => $sff_file,
									    version => $version,
									    version_subdirectory => $sub_dir,
									    );
    ok($read_separation->execute,'execute '. $read_separation->command_name);
}
