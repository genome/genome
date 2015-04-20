package Genome::Config::AnalysisProject::Command::AddConfigBase;

use strict;
use warnings;

use Genome;
use Carp qw();

class Genome::Config::AnalysisProject::Command::AddConfigBase {
    is_abstract => 1,
    is => 'Genome::Config::AnalysisProject::Command::Base',
    has_input => [
        reprocess_existing => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Reprocess any existing instrument data with the new config',
        },
    ],
    has_optional_input => [
        tag => {
            is => 'Genome::Config::Tag',
            is_many => 1,
            doc => 'Tags to associate with the new configuration',
        }
    ]
};

sub valid_statuses {
    return ("Pending", "Hold", "In Progress", "Template");
}

sub execute {
    my $self = shift;

    my @new_items = $self->_create_profile_items();

    $self->_apply_tags(@new_items);

    if($self->reprocess_existing){
        Genome::Config::AnalysisProject::Command::Reprocess->execute(
            analysis_project => $self->analysis_project
        );
    }

    return scalar(@new_items);
}

sub _create_profile_items {
    my $self = shift;
    my $name = $self->command_name;
    Carp::confess sprintf('Command `%s` is missing a way to create profile items.', $name);
}

sub _apply_tags {
    my ($self, @profile_items) = @_;
    for my $profile_item (@profile_items) {
        for my $tag ($self->tag){
            $profile_item->add_tag($tag);
        }
    }

    return 1;
}

1;

