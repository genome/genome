#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT}               = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{NO_LSF} = 1;
}

use Test::More;
use above "Genome";

use Genome::Utility::Test;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ImportedVariationList;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::SoftwareResult::User;

my $TEST_DATA_VERSION = 2;

my $pkg = 'Genome::InstrumentData::Composite::Workflow';
use_ok($pkg) or die('test cannot continue');

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, $TEST_DATA_VERSION);
my $tmp_dir = Genome::Sys->create_temp_directory();
for my $file (qw/indels.hq.vcf/) {
    Genome::Sys->create_symlink(File::Spec->join($data_dir,$file), File::Spec->join($tmp_dir, $file));
}

my $instrument_data_1 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    flow_cell_id => '12345ABXX',
    lane => '1',
    subset_name => '1',
    run_name => 'example',
    id => '-23',
);
my $instrument_data_2 = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    library_id => $instrument_data_1->library_id,
    flow_cell_id => '12345ABXX',
    lane => '2',
    subset_name => '2',
    run_name => 'example',
    id => '-24',
);

my $sample_2 = Genome::Sample->create(
    name => 'sample2',
    id => '-101',
);

my $instrument_data_3 = Genome::InstrumentData::Solexa->__define__(
    flow_cell_id => '12345ABXX',
    lane => '3',
    subset_name => '3',
    run_name => 'example',
    id => '-28',
    sample => $sample_2,
);

my @one_instrument_data = ($instrument_data_3);
my @two_instrument_data = ($instrument_data_1, $instrument_data_2);
my @three_instrument_data = (@two_instrument_data, @one_instrument_data);

my $ref = Genome::Model::Build::ReferenceSequence->get_by_name('GRCh37-lite-build37');

my @alignment_result_one_inst_data = @{construct_alignment_results(\@one_instrument_data, $ref)};
my @alignment_result_two_inst_data = @{construct_alignment_results(\@two_instrument_data, $ref)};
my @alignment_result_three_inst_data = (@alignment_result_one_inst_data, @alignment_result_two_inst_data);

my $merge_result_one_inst_data = construct_merge_result(\@one_instrument_data, $ref);
my $merge_result_two_inst_data = construct_merge_result(\@two_instrument_data, $ref);

my $variation_list_build = construct_variation_list($tmp_dir);

my ($realigner_result_one_inst_data, $recalibrator_result_one_inst_data) = construct_gatk_results(
    $ref, $variation_list_build, $merge_result_one_inst_data
);
my ($realigner_result_two_inst_data, $recalibrator_result_two_inst_data) = construct_gatk_results(
    $ref, $variation_list_build, $merge_result_two_inst_data
);

my $clip_overlap_result_one_inst_data = construct_clip_overlap_result(
    $ref, $merge_result_one_inst_data
);
my $clip_overlap_result_two_inst_data = construct_clip_overlap_result(
    $ref, $merge_result_two_inst_data
);

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build => $ref,
);


subtest 'simple alignments' => sub {
    my $log_directory = Genome::Sys->create_temp_directory();
    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@two_instrument_data,
            ref => $ref,
            force_fragment => 0,
            result_users => $result_users,
        },
        strategy => 'inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::] api v1',
        log_directory => $log_directory,
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments'
    );


    ok($ad->execute, 'executed dispatcher for simple alignments');

    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply([sort @alignment_result_two_inst_data], [sort @ad_results], 'found expected alignment results');
};

subtest 'simple alignments with merge' => sub {
   my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@two_instrument_data,
            ref => $ref,
            force_fragment => 0,
            result_users => $result_users,
        },
        strategy => 'inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::] then merged using picard 1.29 then deduplicated using picard 1.29 api v1',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments with merge'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments with merge');
    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply([sort @alignment_result_two_inst_data, $merge_result_two_inst_data], [sort @ad_results], 'found expected alignment and merge results');
    check_result_bam(@ad_results);
};

subtest "simple alignments of different samples with merge" => sub {
    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@three_instrument_data,
            ref => $ref,
            force_fragment => 0,
            result_users => $result_users,
        },
        strategy => 'inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::] then merged using picard 1.29 then deduplicated using picard 1.29 api v1',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments of different samples with merge'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments of different samples with merge');
    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply([sort @alignment_result_three_inst_data, $merge_result_two_inst_data, $merge_result_one_inst_data], [sort @ad_results], 'found expected alignment and merge results');
    check_result_bam(@ad_results);
};

subtest "simple alignments of different samples with merge and gatk refine" => sub {
    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@three_instrument_data,
            ref => $ref,
            force_fragment => 0,
            variant_list => [$variation_list_build],
            result_users => $result_users,
        },
        strategy => '
            inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::]
            then merged using picard 1.29
            then deduplicated using picard 1.29
            then refined to variant_list using gatk-best-practices 2.4 [-et NO_ET]
            api v1
        ',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments of different samples with merge and gatk refine'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments of different samples with merge and gatk refine');

    my @gatk_results = Genome::InstrumentData::Gatk::BaseRecalibratorBamResult->get(
        bam_source => [$realigner_result_one_inst_data, $realigner_result_two_inst_data],
    );
    is_deeply(
        [sort $recalibrator_result_one_inst_data, $recalibrator_result_two_inst_data],
        [sort @gatk_results],
        'gatk results as expected'
    );

    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply(
        [sort @alignment_result_three_inst_data, @gatk_results],
        [sort @ad_results],
        'found expected alignment and gatk results'
    );
    check_result_bam(@ad_results);
};

subtest "simple alignments of different samples with merge and clip overlap" => sub {


    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@three_instrument_data,
            ref => $ref,
            force_fragment => 0,
            variant_list => [$variation_list_build],
            result_users => $result_users,
        },
        strategy => '
            inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::]
            then merged using picard 1.29
            then deduplicated using picard 1.29
            then refined using clip-overlap 1.0.11
            api v1
        ',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments of different samples with merge and clip_overlap refine'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments of different samples with merge and clip_overlap refine');

    my @clip_overlap_results = Genome::InstrumentData::BamUtil::ClipOverlapResult->get(
        bam_source => [$merge_result_one_inst_data, $merge_result_two_inst_data]
    );
    is_deeply(
        [sort $clip_overlap_result_one_inst_data, $clip_overlap_result_two_inst_data],
        [sort @clip_overlap_results],
        'clip overlap results as expected'
    );

    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply(
        [sort @alignment_result_three_inst_data, @clip_overlap_results],
        [sort @ad_results],
        'found expected alignment and clip_overlap results'
    );
    check_result_bam(@ad_results);
};

subtest "simple alignments of different samples with merge, gatk and clip overlap" => sub {
    my $clip_overlap_result_one_inst_data_2 = construct_clip_overlap_result(
        $ref, $recalibrator_result_one_inst_data
    );
    my $clip_overlap_result_two_inst_data_2 = construct_clip_overlap_result(
        $ref, $recalibrator_result_two_inst_data
    );

    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@three_instrument_data,
            ref => $ref,
            force_fragment => 0,
            variant_list => [$variation_list_build],
            result_users => $result_users,
        },
        strategy => '
            inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::]
            then merged using picard 1.29
            then deduplicated using picard 1.29
            then refined to variant_list using gatk-best-practices 2.4 [-et NO_ET]
            then refined using clip-overlap 1.0.11
            api v1
        ',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments of different samples with merge and clip_overlap refine'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments of different samples with merge and clip_overlap refine');

    my @gatk_results = Genome::InstrumentData::Gatk::BaseRecalibratorBamResult->get(
        bam_source => [$realigner_result_one_inst_data, $realigner_result_two_inst_data],
    );
    is_deeply(
        [sort $recalibrator_result_one_inst_data, $recalibrator_result_two_inst_data],
        [sort @gatk_results],
        'gatk results as expected'
    );

    my @clip_overlap_results = Genome::InstrumentData::BamUtil::ClipOverlapResult->get(
        bam_source => \@gatk_results,
    );
    is_deeply(
        [sort $clip_overlap_result_one_inst_data_2, $clip_overlap_result_two_inst_data_2],
        [sort @clip_overlap_results],
        'clip overlap results as expected'
    );

    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply(
        [sort @alignment_result_three_inst_data, @clip_overlap_results],
        [sort @ad_results],
        'found expected alignment and clip_overlap results'
    );
    check_result_bam(@ad_results);
};

subtest "simple alignments of different samples with merge, clip overlap and gatk" => sub {
    my ($realigner_result_one_inst_data_2, $recalibrator_result_one_inst_data_2) = construct_gatk_results(
        $ref, $variation_list_build, $clip_overlap_result_one_inst_data
    );
    my ($realigner_result_two_inst_data_2, $recalibrator_result_two_inst_data_2) = construct_gatk_results(
        $ref, $variation_list_build, $clip_overlap_result_two_inst_data
    );

    my $ad = Genome::InstrumentData::Composite::Workflow->create(
        inputs => {
            inst => \@three_instrument_data,
            ref => $ref,
            force_fragment => 0,
            variant_list => [$variation_list_build],
            result_users => $result_users,
        },
        strategy => '
            inst aligned to ref using bwa 0.5.9 [-t 4 -q 5::]
            then merged using picard 1.29
            then deduplicated using picard 1.29
            then refined using clip-overlap 1.0.11
            then refined to variant_list using gatk-best-practices 2.4 [-et NO_ET]
            api v1
        ',
    );
    isa_ok(
        $ad,
        'Genome::InstrumentData::Composite::Workflow',
        'created dispatcher for simple alignments of different samples with merge and clip_overlap refine'
    );

    ok($ad->execute, 'executed dispatcher for simple alignments of different samples with merge and clip_overlap refine');

    my @clip_overlap_results = Genome::InstrumentData::BamUtil::ClipOverlapResult->get(
        bam_source => [$merge_result_one_inst_data, $merge_result_two_inst_data],
    );
    is_deeply(
        [sort $clip_overlap_result_one_inst_data, $clip_overlap_result_two_inst_data],
        [sort @clip_overlap_results],
        'clip overlap results as expected'
    );

    my @gatk_results = Genome::InstrumentData::Gatk::BaseRecalibratorBamResult->get(
        bam_source => [$realigner_result_one_inst_data_2, $realigner_result_two_inst_data_2],
    );
    is_deeply(
        [sort $recalibrator_result_one_inst_data_2, $recalibrator_result_two_inst_data_2],
        [sort @gatk_results],
        'gatk results as expected'
    );

    my @ad_result_ids = $ad->_result_ids;
    my @ad_results = Genome::SoftwareResult->get(\@ad_result_ids);
    is_deeply(
        [sort @alignment_result_three_inst_data, @gatk_results],
        [sort @ad_results],
        'found expected alignment and gatk results'
    );
    check_result_bam(@ad_results);
};




sub construct_aligner_index {
    my $reference = shift;

    my $aligner_index = Genome::Model::Build::ReferenceSequence::AlignerIndex->__define__(
        'aligner_version' => '0.5.9',
        'aligner_name' => 'bwa',
        'aligner_params' => '',
        'reference_build' => $reference,
    );

    return $aligner_index;
}

sub construct_alignment_results {
    my $instrument_data_ref = shift;
    my $reference = shift;

    my @alignment_results;
    for my $instrument_data (@{$instrument_data_ref}) {
        my $instrument_data_id = $instrument_data->id;
        my $alignment_result = Genome::InstrumentData::AlignmentResult::Bwa->__define__(
            'reference_build_id' => $reference->id,
            'samtools_version' => 'r599',
            'aligner_params' => '-t 4 -q 5::',
            'aligner_name' => 'bwa',
            'aligner_version' => '0.5.9',
            'picard_version' => '1.29',
            'instrument_data_id' => $instrument_data_id,
            output_dir => '/tmp/fake',
        );
        $alignment_result->lookup_hash($alignment_result->calculate_lookup_hash());
        push @alignment_results, $alignment_result;
    }

    return \@alignment_results;
}

sub construct_variation_list {
    my $directory = shift;

    my $indel_result = Genome::Model::Tools::DetectVariants2::Result::Manual->__define__(
        output_dir => $tmp_dir,
        original_file_path => File::Spec->join($directory, 'indels.hq.vcf'),
    );
    my $variation_list_build = Genome::Model::Build::ImportedVariationList->__define__(
        indel_result => $indel_result,
    );

    return $variation_list_build;
}

sub construct_merge_result {
    my $id_ref = shift;
    my $reference = shift;

    my @id = @{$id_ref};

    my $merge_result = Genome::InstrumentData::AlignmentResult::Merged->__define__(
        aligner_name => 'bwa',
        aligner_version => '0.5.9',
        aligner_params => '-t 4 -q 5::',
        samtools_version => 'r599',
        picard_version => '1.29',
        reference_build_id => $reference->id,
        merger_name => 'picard',
        merger_version => '1.29',
        duplication_handler_name => 'picard',
        duplication_handler_version => '1.29',
        output_dir => '/tmp/fake',
    );
    for my $i (0..$#id) {
        $merge_result->add_input(
            name => 'instrument_data_id-' . $i,
            value_id => $id[$i]->id,
        );
    }
    $merge_result->add_param(
        name => 'instrument_data_id_count',
        value_id=> scalar(@id),
    );
    $merge_result->add_param(
        name => 'instrument_data_id_md5',
        value_id => Genome::Sys->md5sum_data(join(':', sort(map($_->id, @id))))
    );

    $merge_result->add_param(
        name => 'filter_name_count',
        value_id => 0,
    );
    $merge_result->add_param(
        name => 'instrument_data_segment_count',
        value_id => 0,
    );
    $merge_result->lookup_hash($merge_result->calculate_lookup_hash());

    return $merge_result;
}

sub construct_gatk_results {
    my $reference = shift;
    my $variation_list_build = shift;
    my $previous_result = shift;

    my $realigner_result = Genome::InstrumentData::Gatk::IndelRealignerResult->__define__(
        reference_build => $reference,
        bam_source => $previous_result,
        version => 2.4,
        output_dir => '/tmp/fake',
    );
    $realigner_result->add_input(
        name => 'known_sites-1',
        value_obj => $variation_list_build,
    );
    $realigner_result->lookup_hash($realigner_result->calculate_lookup_hash());

    my $recalibrator_result = Genome::InstrumentData::Gatk::BaseRecalibratorBamResult->__define__(
        reference_build => $reference,
        bam_source => $realigner_result,
        version => 2.4,
        output_dir => '/tmp/fake',
    );
    $recalibrator_result->add_input(
        name => 'known_sites-1',
        value_obj => $variation_list_build,
    );
    $recalibrator_result->lookup_hash($recalibrator_result->calculate_lookup_hash());

    return ($realigner_result, $recalibrator_result);
}

sub construct_clip_overlap_result {
    my $reference = shift;
    my $previous_result = shift;

    my $clip_overlap_result = Genome::InstrumentData::BamUtil::ClipOverlapResult->__define__(
        reference_build => $reference,
        bam_source => $previous_result,
        version => '1.0.11',
        output_dir => '/tmp/fake',
    );
    $clip_overlap_result->lookup_hash($clip_overlap_result->calculate_lookup_hash());

    return $clip_overlap_result;
}

sub check_result_bam {
    my @results = @_;

    for my $result (@results) {
        next if ($result->isa("Genome::InstrumentData::AlignmentResult::Bwa"));
        ok($result->bam_file, "Got a bam file from the " . $result->class . " result");
    }
}

done_testing();
