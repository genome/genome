#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';

use File::Spec;
use Test::More tests => 5;

use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Sample;
use Genome::Utility::Test qw(compare_ok);


my $class = 'Genome::InstrumentData::Command::AlignmentResult::Import';
use_ok($class);

my $temp_dir = Genome::Sys->create_temp_directory;

my $fake_bam_path = File::Spec->join($temp_dir, 'fake.bam');
Genome::Sys->write_file($fake_bam_path, "this is a bam file! really!\n");

my $fake_sample = Genome::Test::Factory::Sample->setup_object;
my $fake_reference = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_reference_sequence_build;

my $cmd = $class->create(
    alignment_file => $fake_bam_path,
    sample => $fake_sample,
    reference_build => $fake_reference,
    description => 'import command test',
);
isa_ok($cmd, $class, 'created command');

my $result = $cmd->execute;
isa_ok($result, 'Genome::InstrumentData::AlignmentResult::Merged::External', 'created result');

my $imported_bam = $result->bam_path;
ok(-e $imported_bam, 'file was imported successfully');

compare_ok($imported_bam, $fake_bam_path, 'copied file matches');
