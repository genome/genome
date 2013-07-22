#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

# This test was auto-generated because './Model/SomaticValidation/Command/ValidateSvs/GenerateMergedAssemblies.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::Model::SomaticValidation::Command::ValidateSvs::GenerateMergedAssemblies');

my $temp_dir = Genome::Sys->create_temp_directory();
my $one_line_file = File::Spec->join($temp_dir, 'one_line_file.txt');
Genome::Sys->write_file($one_line_file, "some header content\n");

my $two_line_file = File::Spec->join($temp_dir, 'two_line_file.txt');
Genome::Sys->write_file($two_line_file, "some header content\nsome data and stuff\n");

ok(!Genome::Model::SomaticValidation::Command::ValidateSvs::GenerateMergedAssemblies::_found_somatic_SVs($one_line_file),
        'Got false for a one line file');
ok(Genome::Model::SomaticValidation::Command::ValidateSvs::GenerateMergedAssemblies::_found_somatic_SVs($two_line_file),
        'Got true for a two line file');
