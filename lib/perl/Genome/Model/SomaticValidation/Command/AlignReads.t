#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Test::More tests => 6;

use File::Spec;
use Test::Exception;

use Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::SomaticValidation;


my $class = 'Genome::Model::SomaticValidation::Command::AlignReads';
use_ok($class);

#test shortcutting on prealigned data
my $somval = Genome::Test::Factory::Model::SomaticValidation->setup_object;
my $build = Genome::Model::Build->create(
    model_id => $somval->id,
    data_directory => Genome::Sys->create_temp_directory,
);
Genome::Sys->create_directory(File::Spec->join($build->data_directory, 'alignments'));

my $cmd = $class->create(build_id => $build->id);
isa_ok($cmd, $class, 'created command');

lives_and( sub { ok !$cmd->shortcut }, 'shortcut does not succeed, but does not crash with no premade data');

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object;
my $alignment = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
    output_dir => Genome::Sys->create_temp_directory,
);
$alignment->add_input(
    name => 'instrument_data_id-1',
    value_id => $instrument_data->id,
);
$alignment->add_param(
    name => 'instrument_data_id_count',
    value_id => 1,
);
$alignment->add_param(
    name => 'instrument_data_id_md5',
    value_id => Genome::Sys->md5sum_data($instrument_data->id)
);

$build->add_prealigned_data($alignment);

dies_ok( sub { $cmd->shortcut }, 'shortcut dies with prealigned data for non-matching sample');

$build->tumor_sample($instrument_data->sample);

lives_and( sub { ok $cmd->shortcut }, 'shortcut succeededs with matching prealigned data');

is($build->merged_alignment_result, $alignment, 'prealigned data gets linked as usual');
