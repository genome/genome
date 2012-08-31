package Genome::Model::Command::Define::GenePrediction::Helper;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::Define::GenePrediction::Helper {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    is_abstract => 1,
    has => [
        taxon_id => {
            is => 'Number',
            doc => 'ID of taxon to be used, this can be derived from the assembly model',
        },
        taxon => {
            is => 'Genome::Taxon',
            id_by => 'taxon_id',
            doc => 'Taxon that will be used as the subject of this model',
        },
    ],
    has_optional => [
        assembly_contigs_file => {
            is => 'Path',
            doc => 'Path to the contigs created by an assembly process. If not given, the start_assembly_build, ' .
                   'create_assembly_model, and assembly_processing_profile_name params are used to attempt to ' .
                   'generate an assembly',
        },
        start_assembly_build => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, an assembly build is started if a completed build is not found on the assembly model',
        },
        create_assembly_model => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, an assembly model is created if one cannot be found with the supplied taxon',
        },
        assembly_processing_profile_name => {
            is => 'Number',
            default => 'Velvet Solexa BWA Qual 10 Filter Length 35',
            doc => 'The processing profile used to create assembly models, if necessary',
        },
        subject_name => {
            calculate_from => ['taxon'],
            calculate => q{
                return unless $taxon;
                return $taxon->name;
            },
        },
        subject_type => {
            is => 'Text',
            default => 'species_name',
        },
        subject_class_name => {
            is => 'Text',
            default => 'Genome::Taxon',
        },
    ],
};

sub execute {
    my $self = shift;

    if ($self->taxon and not $self->subject) {
        $self->subject($self->taxon);
    }

    return $self->SUPER::_execute_body();
}

1;

