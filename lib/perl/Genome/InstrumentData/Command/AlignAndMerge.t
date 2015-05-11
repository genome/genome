#!/usr/bin/env genome-perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::More;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::SoftwareResult::User;
use Sub::Override;

my $pkg = 'Genome::InstrumentData::Command::AlignAndMerge';
use_ok($pkg);

my $test_data_dir = __FILE__.'.d';

my $ref_seq_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object;
my $ref_seq_build = Genome::Test::Factory::Build->setup_object(model_id => $ref_seq_model->id);
use Genome::Model::Build::ReferenceSequence;
my $override = Sub::Override->new(
    'Genome::Model::Build::ReferenceSequence::full_consensus_path',
    sub { return File::Spec->join($test_data_dir, 'human_g1k_v37_20_42220611-42542245.fasta'); }
);
use Genome::InstrumentData::AlignmentResult;
my $override2 = Sub::Override->new(
    'Genome::InstrumentData::AlignmentResult::_prepare_reference_sequences',
    sub { return 1; }
);
my $override3 = Sub::Override->new(
    'Genome::InstrumentData::AlignmentResult::filter_non_database_objects',
    sub { my $self = shift; return @_; }
);

my $instrument_data_1 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    flow_cell_id => '12345ABXX',
    lane => '1',
    subset_name => '1',
    run_name => 'example',
    id => '-23',
);
$instrument_data_1->bam_path(File::Spec->join($test_data_dir, '-533e0bb1a99f4fbe9e31cf6e19907133.bam'));
my $instrument_data_2 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    library_id => $instrument_data_1->library_id,
    flow_cell_id => '12345ABXX',
    lane => '2',
    subset_name => '2',
    run_name => 'example',
    id => 'NA12878',
);
$instrument_data_2->bam_path(File::Spec->join($test_data_dir, 'NA12878.20slice.30X.bam'));
my @two_instrument_data = ($instrument_data_1, $instrument_data_2);

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build => $ref_seq_build,
);

my $command = Genome::InstrumentData::Command::AlignAndMerge->create(
    # instrument_data => [@two_instrument_data],
    instrument_data => [$instrument_data_2],
    reference_sequence_build => $ref_seq_build,
    name => 'speedseq',
    version => 'test',
    params => {},
    result_users => $result_users,
    picard_version => '1.46',
    samtools_version => 'r963',
);
ok($command->execute, 'Command executed correctly');
ok($command->result, 'Merged result created');

my $per_lane_result = Genome::InstrumentData::AlignmentResult::Speedseq->get(instrument_data => $command->instrument_data);
ok($per_lane_result, 'Pre lane result created correctly');

done_testing;
