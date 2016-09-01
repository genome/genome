package Genome::Model::ProteinAnnotation;

use strict;
use warnings;
use Genome;

class Genome::Model::ProteinAnnotation {
    is => 'Genome::Model',
    has => [
        subject => {
            is => 'Genome::Taxon',
            id_by => 'subject_id',
            doc => 'taxon of species from which input gene prediction originated',
        },
    ],
    has_param => [
        strategy => {
            is => 'Text',
            doc => 'a string describing what predictors should be run and what should be done with their results',
        },
        chunk_size => {
            is => 'Number',
            doc => 'number of max sequences per fasta chunk should the fasta need to be divided up',
        },
        dump_predictions_to_file => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
        },
        dump_predictions_to_biosql => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
        },
    ],
    has_input => [
        input_fasta_file => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'input_fasta_file', value_class_name => 'UR::Value::Text' ], 
            is_mutable => 1,
        },
    ],
    has_transient => [
        _workflow_inputs => {
            is_optional => 1,
        },
    ],
    doc => 'execute protein prediction and uploads results to a database',
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    unless(-e $self->input_fasta_file && -f $self->input_fasta_file) {
        push @errors,UR::Object::Tag->create(
            type => 'error',
            properties => ['input_fasta_file'],
            desc => 'input_fasta_file does not exist or is not a file'
        );
    }

    return @errors;
}

sub create {
    die(__PACKAGE__ . ' is deprecated.');
}

1;

