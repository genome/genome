#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::SingleSampleGenotype;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq;

use Test::More tests => 7;

my $pkg = 'Genome::Model::Build::SingleSampleGenotype';
use_ok($pkg);

my $test_dir = __FILE__ . '.d';

my $instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object();
my $model = Genome::Test::Factory::Model::SingleSampleGenotype->setup_object();
$model->add_instrument_data($instrument_data);

my $reference_sequence_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object();
my $reference_sequence_build = Genome::Test::Factory::Build->setup_object(model_id => $reference_sequence_model->id);

my $build = $pkg->create(
    model_id => $model->id,
);
isa_ok($build, $pkg, 'created build');

is($build->instrument_data, $instrument_data, 'build has instrument data assigned');

add_merged_alignment_result_to_build($build, 'bam1');

my $other_build = $pkg->create(
    model_id => $model->id,
);
isa_ok($other_build, $pkg, 'created build');
add_merged_alignment_result_to_build($other_build, 'bam2');

my %expected_merged_alignment_result_diff = (
    merged_alignment_result => sprintf(
        'files are not the same (diff -u %s %s)',
        $build->merged_alignment_result->bam_file,
        $other_build->merged_alignment_result->bam_file
    ),
);
is_deeply({$build->_compare_merged_alignment_result($other_build)}, \%expected_merged_alignment_result_diff, 'Merged alignment result diffs found');

sub add_merged_alignment_result_to_build {
    my $build = shift;
    my $id = shift;

    my $speedseq_result = Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq->setup_object(
        instrument_data => [$build->instrument_data],
        reference_build => $reference_sequence_build,
        id => $id,
        output_dir => $test_dir,
        map {; $_ => $build->model->$_ } (qw(aligner_name aligner_version aligner_params)),
    );
    $speedseq_result->add_user(user => $build, label => 'merged_alignment_result');

    is($build->merged_alignment_result, $speedseq_result, 'alignment result added to build');
}
