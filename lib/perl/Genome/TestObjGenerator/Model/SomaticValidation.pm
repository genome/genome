package Genome::TestObjGenerator::Model::SomaticValidation;
use Genome::TestObjGenerator::Model;
@ISA = (Genome::TestObjGenerator::Model);

use strict;
use warnings;
use Genome;
use Genome::TestObjGenerator::ProcessingProfile::SomaticValidation;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::SomaticValidation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::TestObjGenerator::ProcessingProfile::SomaticValidation->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Individual->create(name => Genome::TestObjGenerator::Util::generate_name("test_individual"));
}

# This isn't possible to due within the constraints of the current
# TestObjGenerator API so at least pulling this up here for better
# reuse/refactoring.
my $count;
sub setup_somatic_variation_build {
    use Genome::TestObjGenerator::Build;
    use Genome::TestObjGenerator::Individual;
    use Genome::TestObjGenerator::Model::ReferenceAlignment;
    use Genome::TestObjGenerator::Model::SomaticVariation;
    use Genome::TestObjGenerator::Sample;
    use Genome::TestObjGenerator::InstrumentData::Solexa;
    use Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment;
    use Genome::TestObjGenerator::ProcessingProfile::SomaticVariation;

    $count++;

    my $test_profile = Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment->setup_object(
        sequencing_platform => 'solexa',
        dna_type => 'cdna',
        read_aligner_name => 'bwa',
        snv_detection_strategy => 'samtools [--test ' . $count . ']', # to make each processing profile unique
    );

    my $test_individual = Genome::TestObjGenerator::Individual->setup_object();

    my $test_sample = Genome::TestObjGenerator::Sample->setup_object(
        source_id => $test_individual->id,
    );

    my $test_control_sample = Genome::TestObjGenerator::Sample->setup_object(
        source_id => $test_individual->id,
    );

    my $test_instrument_data = Genome::TestObjGenerator::InstrumentData::Solexa->setup_object();

    my $test_model = Genome::TestObjGenerator::Model::ReferenceAlignment->setup_object(
        subject_name => $test_sample->name,
        processing_profile_id => $test_profile->id,
    );

    my $add_ok = $test_model->add_instrument_data($test_instrument_data);

    my $test_build = Genome::TestObjGenerator::Build->setup_object(
        model_id => $test_model->id,
    );

    my $test_model_two = Genome::TestObjGenerator::Model::ReferenceAlignment->setup_object(
        subject_name => $test_control_sample->name,
        processing_profile_id => $test_profile->id,
    );

    $add_ok = $test_model_two->add_instrument_data($test_instrument_data);

    my $test_build_two = Genome::TestObjGenerator::Build->setup_object(
        model_id => $test_model_two->id,
    );

    my $test_somvar_pp = Genome::TestObjGenerator::ProcessingProfile::SomaticVariation->setup_object(
        snv_detection_strategy => 'samtools r599 [--test ' . $count . ']',
        tiering_version => 1,
    );

    my $somvar_model = Genome::TestObjGenerator::Model::SomaticVariation->setup_object(
        tumor_model => $test_model,
        normal_model => $test_model_two,
        processing_profile_id => $test_somvar_pp->id,
    );

    my $somvar_build = Genome::TestObjGenerator::Build->setup_object(
        tumor_build => $test_build_two,
        normal_build => $test_build,
        model_id => $somvar_model->id,
    );

    return $somvar_build;
}

1;
