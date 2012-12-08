package Genome::DruggableGene::Command::DrugNameGroup::RemoveAll;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::DrugNameGroup::RemoveAll {
    is => 'Genome::Command::Base',
    doc => 'Remove all DrugNameGroups (for the purpose of rebuilding them after a new import)'
};

sub help_brief {'Remove druggable gene DrugNameGroups for the purpose of rebuilding them'}

sub help_synopsis { help_brief }

sub help_detail { help_brief }

sub execute {
    my $self = shift;
    my @bridges = Genome::DruggableGene::DrugNameGroupBridge->get();
    my @groups = Genome::DruggableGene::DrugNameGroup->get();
    $self->status_message("Found " . scalar(@groups) . " drug name group entries ... deleting these now");
    map($_->delete, (@bridges, @groups));
    return 1;
}

1;
