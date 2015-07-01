package Genome::Model::Command::Define::GenotypeMicroarray;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Model::Command::Define::GenotypeMicroarray {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has => [
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'Variation list [DBSNP] build for this model.',
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
genome model define genotype-microarray 
  --subject-name SAMPLE_NAME 
  --processing-profile-name PROCESSING_PROFILE_NAME 
  --dbsnp-build GRCh37-lite-build37
EOS
}

sub help_detail {
    return "Define a new genome model with genotype information based on microarray data."
}

sub type_specific_parameters_for_create {
    my $self = shift;
    return ($self->SUPER::type_specific_parameters_for_create, dbsnp_build => $self->dbsnp_build);
}

sub execute {
    my $self = shift;

    my $processing_profile = $self->processing_profile;
    unless ($processing_profile) {
        Carp::confess "Could not resolve processing profile!";
    }

    return $self->SUPER::_execute_body(@_);
}

1;

