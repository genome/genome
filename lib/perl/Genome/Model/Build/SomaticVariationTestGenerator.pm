package Genome::Model::Build::SomaticVariationTestGenerator;

use strict;
use warnings;
use Genome;

sub create_somatic_variation_model {
    my $ref_align_pp = Genome::ProcessingProfile::ReferenceAlignment->create(
        name => 'ref_align_pp',
        sequencing_platform => 'solexa',
        dna_type => 'cdna',
        read_aligner_name => 'bwa',
        snv_detection_strategy => 'samtools',
    );
    my $reference_sequence_build = Genome::Model::Build::ReferenceSequence->get_by_name('NCBI-human-build36');

    create_test_subjects();
    my $tumor_model = Genome::Model->create(
        name => 'test_reference_aligment_model_TUMOR',
        subject_name => 'test_tumor_subject',
        subject_type => 'sample_name',
        processing_profile_id => $ref_align_pp->id,
        reference_sequence_build => $reference_sequence_build,
    );
    my $test_instrument_data = Genome::InstrumentData::Solexa->create();
    $tumor_model->add_instrument_data($test_instrument_data);
    my $tumor_build = Genome::Model::Build->create(
        model_id => $tumor_model->id,
        data_directory => Genome::Sys->create_temp_directory(),
    );
    Genome::Model::Event::Build->__define__(
        build_id => $tumor_build->id,
        event_type => 'genome model build',
        event_status => 'Succeeded',
        model_id => $tumor_model->id,
        date_completed => '1999-01-01 15:19:01',
    );

    my $normal_model = Genome::Model->create(
        name => 'test_reference_aligment_model_mock_NORMAL',
        subject_name => 'test_normal_subject',
        subject_type => 'sample_name',
        processing_profile_id => $ref_align_pp->id,
        reference_sequence_build => $reference_sequence_build,
    );
    $normal_model->add_instrument_data($test_instrument_data);
    my $normal_build = Genome::Model::Build->create(
        model_id => $normal_model->id,
        data_directory => Genome::Sys->create_temp_directory(),
    );
    Genome::Model::Event::Build->__define__(
        build_id => $normal_build->id,
        event_type => 'genome model build',
        event_status => 'Succeeded',
        model_id => $normal_model->id,
        date_completed => '1999-01-01 15:19:01',
    );

    my $annotation_build = Genome::Model::Build::ImportedAnnotation->__define__(model_id => '-1');
    my $somvar_pp = Genome::ProcessingProfile::SomaticVariation->create(
        name => 'somvar_pp',
        snv_detection_strategy => 'samtools r599 [--test=1]',
        tiering_version => 1,
    );
    my $somvar_model = Genome::Model::SomaticVariation->create(
        processing_profile => $somvar_pp,
        tumor_model => $tumor_model,
        normal_model => $normal_model,
        name => 'test somvar model',
        annotation_build => $annotation_build,
        previously_discovered_variations =>
                create_pdv_build($reference_sequence_build),
    );
    return $somvar_model, $tumor_build, $normal_build;
}


# create Samples so they can be used in tumor/normal_models.
sub create_test_subjects {
    # don't create them more than once.
    unless (Genome::Individual->get(name => 'test_individual')) {
        my $test_individual = Genome::Individual->create(
            common_name => 'TEST',
            name => 'test_individual',
        );
        Genome::Sample->create(
            name => 'test_tumor_subject',
            source_id => $test_individual->id,
        );
        Genome::Sample->create(
            name => 'test_normal_subject',
            source_id => $test_individual->id,
        );
    }
}

sub create_pdv_model {
    my ($reference_sequence_build) = @_;

    my $pdv_model = Genome::Model::ImportedVariationList->get(
            name => "test imported-variation-list model",
    );
    unless ($pdv_model) {
        my $pdv_pp = Genome::ProcessingProfile::ImportedVariationList->create(name => 'pdv_pp');
        $pdv_model = Genome::Model::ImportedVariationList->create(
            processing_profile => $pdv_pp,
            subject_name => 'test_normal_subject',
            subject_type => 'sample_name',
            reference => $reference_sequence_build,
            name => "test imported-variation-list model",
        );
    }
    return $pdv_model;
}

sub create_pdv_build {
    my ($reference_sequence_build) = @_;

    create_test_subjects();
    my $pdv_model = create_pdv_model();
    my $pdv_build = Genome::Model::Build::ImportedVariationList->create(
        model_id => $pdv_model->id,
        data_directory => Genome::Sys->create_temp_directory(),
    );
    return $pdv_build;
}


sub setup_test_build {
    my ($somvar_model, $tumor_build, $normal_build) = create_somatic_variation_model();

    my $somvar_build = Genome::Model::Build::SomaticVariation->create(
        model_id => $somvar_model->id,
        data_directory => Genome::Sys->create_temp_directory(),
    );
    return $somvar_build, $somvar_model;
}

sub setup_test_model {
    my ($somvar_model, $tumor_build, $normal_build) = create_somatic_variation_model();

    my $somvar_build = Genome::Model::Build::SomaticVariation->__define__(
        model_id => $somvar_model->id,
        data_directory => Genome::Sys->create_temp_directory(),
        tumor_build => $tumor_build,
        normal_build => $normal_build,
    );
    my $e = Genome::Model::Event::Build->__define__(
        build_id => $somvar_build->id,
        event_type => 'genome model build',
        event_status => 'Succeeded',
        model_id => $somvar_model->id,
        date_completed => '1999-01-01 15:19:01',
    );


    my $output_dir = Genome::Sys->create_temp_directory();
    my $tier_result = Genome::Model::Tools::DetectVariants2::Classify::Tier->__define__(
        output_dir => $output_dir,
        variant_type => 'snv',
    ); 
    $tier_result->lookup_hash($tier_result->calculate_lookup_hash());
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

