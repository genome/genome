#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
use File::Temp 'tempdir';
use Test::More tests => 5;

BEGIN {
        use_ok('Genome::Model::Tools::Hgmi::Predict');
        use_ok('Genome::Model::Tools::Hgmi::DirBuilder');
}

my $test_dir = tempdir(
    'Genome-Model-Tools-Hgmi-Predict-XXXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
    UNLINK => 1,
);

my $d = Genome::Model::Tools::Hgmi::DirBuilder->create(
    path => $test_dir,
    org_dirname => "B_catenulatum",
    assembly_version_name => "Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb",
    assembly_version => "Version_1.0",
    pipe_version => "Version_1.0",
    cell_type => "BACTERIA"
);

isa_ok($d,'Genome::Model::Tools::Hgmi::DirBuilder');
ok($d->execute(), "executed dir builder command");

my $p = Genome::Model::Tools::Hgmi::Predict->create(
    organism_name => "Bifidobacterium_catenulatum",
    locus_tag => "BIFCATDFT",
    project_type => "HGMI",
    work_directory => $test_dir,
    dev => 1
);

isa_ok($p, 'Genome::Model::Tools::Hgmi::Predict');
