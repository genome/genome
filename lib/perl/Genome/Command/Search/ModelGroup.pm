package Genome::Command::Search::ModelGroup;

use strict;
use warnings;

class Genome::Command::Search::ModelGroup {
    is => 'Genome::Command::Search::Base',
};

sub display_single {
    my ($self, $object) = @_;
    my @models = $object->models;
    if (length(@models)) {
        Genome::Model::Command::Status->execute(models => \@models);
    } else {
        warn sprintf("No models found associated with model group %s (%s).\n",
            $object->name, $object->id);
        Genome::ModelGroup::Command::list->execute(model_groups => [$object]);
    }
}

sub display_many {
    my ($self, $objects) = @_;
    die $self->error_message("Cannot display multiple model groups");
}
