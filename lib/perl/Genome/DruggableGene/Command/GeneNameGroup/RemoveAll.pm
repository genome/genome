package Genome::DruggableGene::Command::GeneNameGroup::RemoveAll;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::GeneNameGroup::RemoveAll {
    is => 'Genome::Command::Base',
    doc => 'Remove all GeneNameGroups (for the purpose of rebuilding them after a new entrez import)'
};

sub help_brief {'Remove druggable gene GeneNameGroups for the purpose of rebuilding them'}

sub help_synopsis { help_brief }

sub help_detail { help_brief }

sub execute {
    my $self = shift;
    my @bridges = Genome::DruggableGene::GeneNameGroupBridge->get();
    my @groups = Genome::DruggableGene::GeneNameGroup->get();
    map($_->delete, (@bridges, @groups));
}

1;
