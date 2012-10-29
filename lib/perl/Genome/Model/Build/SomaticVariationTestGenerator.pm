package Genome::Model::Build::SomaticVariationTestGenerator;

use strict;
use warnings;
use Genome;

sub setup_test_model {
    my $test_profile = Genome::ProcessingProfile::ReferenceAlignment->create(
        name => 'test_profile',
        sequencing_platform => 'solexa',
        dna_type => 'cdna',
        read_aligner_name => 'bwa',
        snv_detection_strategy => 'samtools',
    );

    my $test_somvar_pp = Genome::ProcessingProfile::SomaticVariation->create(
        name => 'test somvar pp',
        snv_detection_strategy => 'samtools r599 [--test=1]',
        tiering_version => 1,
    );

    my $test_pdv_profile = Genome::ProcessingProfile::ImportedVariationList->create(
        name => 'test_pdv_profile',
    );

    my $annotation_build = Genome::Model::Build::ImportedAnnotation->__define__(
        model_id => '-1',
    );

    my $test_individual = Genome::Individual->create(
        common_name => 'TEST',
        name => 'test_individual',
    );

    my $test_sample = Genome::Sample->create(
        name => 'test_subject',
        source_id => $test_individual->id,
    );

    my $test_control_sample = Genome::Sample->create(
        name => 'test_control_subject',
        source_id => $test_individual->id,
    );

    my $test_instrument_data = Genome::InstrumentData::Solexa->create(
    );

    my $reference_sequence_build = Genome::Model::Build::ReferenceSequence->get_by_name('NCBI-human-build36');

    my $test_model = Genome::Model->create(
        name => 'test_reference_aligment_model_TUMOR',
        subject_name => 'test_subject',
        subject_type => 'sample_name',
        processing_profile_id => $test_profile->id,
        reference_sequence_build => $reference_sequence_build,
    );

    my $add_ok = $test_model->add_instrument_data($test_instrument_data);

    my $temp_build_data_dir = Genome::Sys->create_temp_directory;

    my $test_build = Genome::Model::Build->create(
        model_id => $test_model->id,
        data_directory => $temp_build_data_dir,
    );

    my $test_model_two = Genome::Model->create(
        name => 'test_reference_aligment_model_mock_NORMAL',
        subject_name => 'test_control_subject',
        subject_type => 'sample_name',
        processing_profile_id => $test_profile->id,
        reference_sequence_build => $reference_sequence_build,
    );

    $add_ok = $test_model_two->add_instrument_data($test_instrument_data);

    my $test_build_two = Genome::Model::Build->create(
        model_id => $test_model_two->id,
        data_directory => $temp_build_data_dir,
    );

    my $pdv_model = Genome::Model::ImportedVariationList->create(
        processing_profile => $test_pdv_profile,
        subject_name => 'test_control_subject',
        subject_type => 'sample_name',
        reference => $reference_sequence_build,
        name => "test imported-variation-list model",
    );

    my $temp_pdv_build_data_dir1 = Genome::Sys->create_temp_directory;
    my $pdv_build_1 = Genome::Model::Build::ImportedVariationList->create(
        model_id => $pdv_model->id,
        data_directory => $temp_pdv_build_data_dir1,
    );

    my $somvar_model = Genome::Model::SomaticVariation->create(
        tumor_model => $test_model,
        normal_model => $test_model_two,
        name => 'test somvar model',
        processing_profile => $test_somvar_pp,
        annotation_build => $annotation_build,
        previously_discovered_variations => $pdv_build_1,
    );

    my $somvar_build = Genome::Model::Build::SomaticVariation->__define__(
        model_id => $somvar_model->id,
        data_directory => $temp_build_data_dir,
        tumor_build => $test_build_two,
        normal_build => $test_build,
    );
    my $e = Genome::Model::Event::Build->__define__(
        build_id => $somvar_build->id,
        event_type => 'genome model build',
        event_status => 'Succeeded',
        model_id => $somvar_model->id,
        date_completed => '1999-01-01 15:19:01',
    );


    my $output_dir = Genome::Sys->create_temp_directory;
    my $tier_result = Genome::Model::Tools::DetectVariants2::Classify::Tier->__define__(
        output_dir => $output_dir,
        variant_type => 'snv',
    );
    $tier_result->add_user(user => $somvar_build, label => 'uses');

    #no relation to any real tiers
    my $tier1_data = <<EOBED
1	23	24	A/T
1	29	30	A/T
2	320	321	A/T
EOBED
;
    my $tier2_data = <<EOBED
1	950	951	A/T
1	1005	1006	A/T
EOBED
;
    my $tier3_data = <<EOBED
2	98	99	A/T
2	101	102	A/T
2	500	501	A/T
EOBED
;
    my $tier4_data = <<EOBED
3	1	2	A/T
EOBED
;

    my $tier1_file = join('/', $output_dir, 'tier1.hq.bed');
    my $tier2_file = join('/', $output_dir, 'tier2.hq.bed');
    my $tier3_file = join('/', $output_dir, 'tier3.hq.bed');
    my $tier4_file = join('/', $output_dir, 'tier4.hq.bed');

    Genome::Sys->write_file($tier1_file, $tier1_data);
    Genome::Sys->write_file($tier2_file, $tier2_data);
    Genome::Sys->write_file($tier3_file, $tier3_data);
    Genome::Sys->write_file($tier4_file, $tier4_data);

    return $somvar_model;
}

1;

