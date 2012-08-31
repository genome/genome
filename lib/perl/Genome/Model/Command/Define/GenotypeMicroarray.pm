package Genome::Model::Command::Define::GenotypeMicroarray;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Model::Command::Define::GenotypeMicroarray {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            id_by => 'reference_id',
            doc => 'reference sequence build for this model',
        },
        reference_id => {
            is => 'Text',
            is_input => 1,
            doc => 'id of reference build',
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
genome model define genotype-microarray 
  --subject-name SAMPLE_NAME 
  --processing-profile-name PROCESSING_PROFILE_NAME 
  --reference GRCh37-lite-build37
EOS
}

sub help_detail {
    return "Define a new genome model with genotype information based on microarray data."
}

sub type_specific_parameters_for_create {
    my $self = shift;
    return (reference_sequence_build => $self->reference);
}

sub execute {
    my $self = shift;

    my $processing_profile = $self->processing_profile;
    unless ($processing_profile) {
        Carp::confess "Could not resolve processing profile!";
    }

    unless ($self->processing_profile->name =~ /wugc/) {
        Carp:confess "GenotypeMicroarray Models must use one of the [microarray-type]/wugc processing-profiles.";
    }

    return $self->SUPER::_execute_body(@_);
}

1;

