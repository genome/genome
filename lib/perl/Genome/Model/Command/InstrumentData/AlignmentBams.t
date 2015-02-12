#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::DiskAllocation;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ReferenceAlignment;

use Test::More tests => 5;

use constant NUM_INSTRUMENT_DATA => 3;

my $pkg = 'Genome::Model::Command::InstrumentData::AlignmentBams';
use_ok($pkg);

my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();

for my $i (1..NUM_INSTRUMENT_DATA) {
    my $inst_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
        lane => $i,
        flow_cell_id => 'TEST1AAXX',
    );
    $model->add_instrument_data($inst_data);
}


my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);
for my $input ($build->instrument_data_inputs) {

    my ($params) = $model->processing_profile->params_for_alignment($input);

    my $ar = Genome::InstrumentData::AlignmentResult::Bwa->__define__(
        %$params,
        output_dir => '/fake/' . $input->value_id,

    );
    $ar->lookup_hash($ar->calculate_lookup_hash);
    $ar->add_user(label => 'created', user => $build);

    my $qc = Genome::InstrumentData::AlignmentResult::Merged::BamQc->__define__(
        alignment_result_id => $ar->id,
        output_dir => '/qc-fake/' . $input->value_id,
    );
    $qc->add_user(label => 'created', user => $build);
    Genome::Test::Factory::DiskAllocation->setup_object(owner => $qc, creation_time => '1970-01-01T00:00:00Z');
}

my $dir = Genome::Sys->create_temp_directory();
my $cmd = $pkg->create(
    build => $build,
    outdir => $dir,
);
isa_ok($cmd, $pkg);

ok($cmd->execute, 'executed command');

my @files = glob(File::Spec->join($dir, '*'));
is(scalar(@files), 1, 'wrote output file') or die 'cannot continue without output';
my @data = Genome::Sys->read_file($files[0]);

is(scalar(@data), NUM_INSTRUMENT_DATA+1, 'one line per instrument data plus a header');
