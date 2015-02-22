package Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::Model;
@ISA = (Genome::Test::Factory::Model);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::Build;
use Genome::Test::Factory::ProcessingProfile::ReferenceSequence;

our @required_params = qw(subject);

sub generate_obj {
    my $self = shift;

    my $m = Genome::Model::ReferenceSequence->create(@_);
    return $m;
}

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::ReferenceSequence->setup_object();
    return $p->id;
}

sub create_subject {
    return Genome::Sample->create(name => Genome::Test::Factory::Util::generate_name("test_subject"));
}

sub setup_reference_sequence_build {
    my $self = shift;
    my $fasta_file = shift;

    my $test_model = Genome::Test::Factory::Model::ReferenceSequence->setup_object(
        fasta_file => $fasta_file,
    );

    my $test_build = Genome::Test::Factory::Build->setup_object(
        model_id => $test_model->id,
        status => 'Succeeded',
    );
    return $test_build;
}

1;
