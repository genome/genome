package Genome::Model::ReferenceSequence::Command::Converter::Create;

use strict;
use warnings;

use Genome;


class Genome::Model::ReferenceSequence::Command::Converter::Create {
    is => 'Command::V2',
    has_input => [
        source_reference => { is => 'Genome::Model::Build::ReferenceSequence', doc => 'the reference to convert from' },
        destination_reference => { is => 'Genome::Model::Build::ReferenceSequence', doc => 'the reference to convert to' },
        algorithm => { is => 'Text', doc => 'method to do the conversion (valid values are methods of Genome::FeatureList::Converter)' },
    ],
    has_optional_input => [
        resource => { is => 'Text', doc => 'additional resource for the algorithm to use (e.g. chain file for liftOver)' },
    ],
};

sub help_brief {
    "Create a new reference-sequence converter.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt model reference-sequence converter create --source NCBI-human-build36 --destination NCBI-human-build36 --algorithm no_op
EOS
}

sub help_detail {
    return <<EOS
Create a new reference-sequence converter.  This specifies a way to convert coordinates in one reference to coordinates in another.
EOS
}

sub execute {
    my $self = shift;

    my $existing_converter = Genome::Model::Build::ReferenceSequence::Converter->exists_for_references($self->source_reference, $self->destination_reference);
    if($existing_converter) {
        $self->error_message('There is already an existing reference-sequence converter for the specified references.');
        return;
    }

    my $converter = Genome::Model::Build::ReferenceSequence::Converter->create(
        source_reference_build => $self->source_reference,
        destination_reference_build => $self->destination_reference,
        algorithm => $self->algorithm,
        ($self->resource ? (resource => $self->resource) : () ),
    );
    unless($converter) {
        $self->error_message('Failed to create reference-sequence converter.');
        return;
    }

    if($converter->__errors__) {
        $self->error_message('Errors creating reference-sequence converter: ' . join(' ',map($_->__display_name__, $self->__errors__)));
        $converter->delete;
        return;
    }

    return 1;
}

1;
