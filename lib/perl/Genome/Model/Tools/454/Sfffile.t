#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

BEGIN {
    my $archos = `uname -a`;
    if ($archos !~ /64/) {
        plan skip_all => "Must run from 64-bit machine";
    }
    plan tests => 5;
    use_ok('Genome::Model::Tools::454::Sfffile');
}

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_file = $tmp_dir .'/tmp.sff';

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-454-Newbler/R_2008_09_22_14_01_00_FLX12345678_TEST_12345678';

my @sff_files = glob($data_dir.'/*.sff');

my @version_subdirs = qw/ offInstrumentApps mapasm454_source /;

foreach my $sub_dir (@version_subdirs) {
    my $version;
    $version = '2.0.00.20-1' if $sub_dir eq 'offInstrumentApps';
    $version = '10282008' if $sub_dir eq 'mapasm454_source';

    my $sfffile = Genome::Model::Tools::454::Sfffile->create(
							     in_sff_files => \@sff_files,
							     out_sff_file => $out_file,
							     version => $version,
							     version_subdirectory => $sub_dir,
							     );
    isa_ok($sfffile,'Genome::Model::Tools::454::Sfffile');
    ok($sfffile->execute,'execute '. $sfffile->command_name);
}
