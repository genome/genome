package Genome::FeatureList::Command::Update;

use strict;
use warnings;

use Genome;

class Genome::FeatureList::Command::Update{
    is => 'Command::V2',
    has_input => [
        feature_list => { 
            is => 'Genome::FeatureList', 
            doc => 'The feature list to update', 
            shell_args_position => 1,
            is_many => 1,
        },
    ],
    has_optional_input => [
        format => {
            is => 'Text',
            doc => 'New format for the feature list',
            valid_values => Genome::FeatureList->__meta__->property('format')->valid_values
        },
        source => {
            is => 'Text',
            doc => 'New source for the feature list',
        },
        reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => "New reference sequence for the feature list",
        },
        subject => {
            is => 'Genome::Model::Build',
            doc => "New subject for the feature list",
        },
        content_type => {
            is => 'Text',
            doc => 'New content type for the feature list',
            valid_values => Genome::FeatureList->__meta__->property('content_type')->valid_values,
        },
        description => {
            is => 'Text',
            doc => 'New description for the feature list',
        },
    ],
};

sub help_brief {

}

sub help_synopsis {

}

sub help_detail{

}

sub execute{
    my $self = shift;
    my @feature_lists = $self->feature_list;
    my @properties = $self->__meta__->properties;
    my @property_names = map { $_->property_name } grep { $_->can('is_input') and $_->is_input } @properties;

    for my $feature_list (@feature_lists) {
        for my $property (@property_names) {
            next if $property eq 'feature_list';
            if ($self->$property){
                my $new_value = $self->$property;
                $self->status_message("setting $property to $new_value\n");
                $feature_list->$property($new_value);
            }
        }
    }

    return 1;
}

1;
