package Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::ProcessingProfile::SomaticValidation;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::SomaticValidation->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::SomaticValidation->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Individual->create(name => Genome::Test::Factory::Util::generate_name("test_individual"));
}

# This isn't possible to due within the constraints of the current
# Test::Factory API so at least pulling this up here for better
# reuse/refactoring.
my $count;
sub setup_somatic_variation_build {
    use Genome::Test::Factory::Build;
    use Genome::Test::Factory::Individual;
    use Genome::Test::Factory::Model::ReferenceAlignment;
    use Genome::Test::Factory::Model::SomaticVariation;
    use Genome::Test::Factory::Sample;
    use Genome::Test::Factory::InstrumentData::Solexa;
    use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment;
    use Genome::Test::Factory::ProcessingProfile::SomaticVariation;

    $count++;

    my $test_profile = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object(
        sequencing_platform => 'solexa',
        dna_type => 'cdna',
        read_aligner_name => 'bwa',
        snv_detection_strategy => 'samtools [--test ' . $count . ']', # to make each processing profile unique
    );

    my $test_individual = Genome::Test::Factory::Individual->setup_object();

    my $test_sample = Genome::Test::Factory::Sample->setup_object(
        source_id => $test_individual->id,
    );

    my $test_control_sample = Genome::Test::Factory::Sample->setup_object(
        source_id => $test_individual->id,
    );

    my $test_instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object();

    my $test_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        subject_id => $test_sample->id,
        processing_profile_id => $test_profile->id,
    );

    my $add_ok = $test_model->add_instrument_data($test_instrument_data);

    my $test_build = Genome::Test::Factory::Build->setup_object(
        model_id => $test_model->id,
        status => 'Succeeded',
    );

    my $test_model_two = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        subject_id => $test_control_sample->id,
        processing_profile_id => $test_profile->id,
    );

    $add_ok = $test_model_two->add_instrument_data($test_instrument_data);

    my $test_build_two = Genome::Test::Factory::Build->setup_object(
        model_id => $test_model_two->id,
        status => 'Succeeded',
    );

    my $test_somvar_pp = Genome::Test::Factory::ProcessingProfile::SomaticVariation->setup_object(
        snv_detection_strategy => 'samtools r599 [--test ' . $count . ']',
        tiering_version => 1,
    );

    my $somvar_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
        tumor_model => $test_model,
        normal_model => $test_model_two,
        processing_profile_id => $test_somvar_pp->id,
    );

    my $somvar_build = Genome::Test::Factory::Build->setup_object(
        model_id => $somvar_model->id,
        status => 'Succeeded',
    );

    return $somvar_build;
}

1;
