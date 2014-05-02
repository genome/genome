package Genome::Annotation::Expert::Vep::Adaptor;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Annotation::Expert::Vep::Adaptor {
    is => "Genome::Annotation::AdaptorBase",

    has_planned_output => [
        ensembl_version => {
            is => 'String',
        },
        species => { is => 'Text', },
        polyphen => { is => 'String', },
        sift => { is => 'String', },
        terms => {is => 'String',},
        regulatory => {is => 'Boolean',},
        canonical => {is => 'Boolean',},
        plugins => {is => 'String',
                    is_many => 1,
                    is_optional => 1},
        feature_list_names_and_tags => {is => 'HASH',},
        plugins_version => { is => 'String', },
    ],
    has_output => [
        feature_list_ids_and_tags => {
            is => 'string',
            is_many => 1,
            is_optional => 1,
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
    ],
};

sub name {
    'vep';
}

sub resolve_expert_specific_attributes_from_build {
    my $self = shift;

    my @strings;
    for my $name (keys %{$self->feature_list_names_and_tags}) {
        my $feature_list = $self->build->get_feature_list($name);
        push @strings, sprintf('%s:%s', $feature_list->id,
            $self->feature_list_names_and_tags->{$name});
    }
    $self->feature_list_ids_and_tags(\@strings);

    $self->reference_build($self->build->reference_sequence_build);
}

1;
