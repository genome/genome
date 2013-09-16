package Genome::Model::Build::SomaticVariationTestGenerator;

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ImportedVariationList;
use Genome::Test::Factory::Build;

sub create_somatic_variation_model {
    
    my $tumor_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(subject_name => "test_subject");
    my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(processing_profile_id => $tumor_model->processing_profile->id, subject_name => $tumor_model->subject_name);
    my $tumor_build = Genome::Test::Factory::Build->setup_object(
        model_id => $tumor_model->id, status => "Succeeded");

    my $reference_sequence_build = $normal_model->reference_sequence_build;
    my $normal_build = Genome::Test::Factory::Build->setup_object(
        model_id => $normal_model->id,
        status => "Succeeded",
    );

    my $annotation_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object(name => "test_annotation_build");
    my $annotation_build = Genome::Test::Factory::Build->setup_object(model_id => $annotation_model->id, version => "1", name => "test_annotation_build/1", status => "Succeeded");

    my $somvar_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(normal_model => $normal_model, tumor_model => $tumor_model, annotation_build => $annotation_build, subject_name => $normal_model->subject_name, previously_discovered_variations => create_pdv_build());
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

    my $pdv_model = Genome::Test::Factory::Model::ImportedVariationList->setup_object(name => "test imported-variation-list model");
    return $pdv_model;
}

sub create_pdv_build {
    my ($reference_sequence_build) = @_;

    my $pdv_model = create_pdv_model();
    my $pdv_build = Genome::Test::Factory::Build->setup_object(
        model_id => $pdv_model->id,
    );
    return $pdv_build;
}

sub setup_test_build {
    my %args = @_;
    my $data_dir;
    my $annot_data_dir;
    my $normal;
    my $tumor;

    if ($args{som_var_dir}) {
        $data_dir = delete $args{som_var_dir};
    }
    else {
        $data_dir = Genome::Sys->create_temp_directory();
    }
    if ($args{annot_dir}) {
        $annot_data_dir = delete $args{annot_dir};
    }
    else {
        $annot_data_dir = Genome::Sys->create_temp_directory();
    }
    
    my ($somvar_model, $tumor_build, $normal_build) = create_somatic_variation_model();

    my $somvar_build = Genome::Test::Factory::Build->setup_object(
        model_id => $somvar_model->id,
        data_directory => $data_dir,
    );
    $somvar_build->annotation_build->data_directory($annot_data_dir);
    return $somvar_build, $somvar_model;
}

sub setup_test_model {
    my ($somvar_model, $tumor_build, $normal_build) = create_somatic_variation_model();

    my $somvar_build = Genome::Test::Factory::Build->setup_object(
        model_id => $somvar_model->id,
        tumor_build => $tumor_build,
        normal_build => $normal_build,
        status => "Succeeded",
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

