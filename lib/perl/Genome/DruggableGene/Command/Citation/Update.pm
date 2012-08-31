package Genome::DruggableGene::Command::Citation::Update;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::Citation::Update {
    is => 'Genome::Command::Base',
    has_input => [
        citation => {
            is => 'Genome::DruggableGene::Citation',
            doc => 'Citation to update',
            shell_args_position => 1,
            is_many => 1,
        },
    ],
    has_optional_input => [
        base_url => {
            is => 'Text',
            doc => 'New base_url for the citation',
        },
    ],
};


sub help_brief {

}

sub help_synopsis {

}

sub help_detail {

}

sub execute {
    my $self = shift;
    my @citations = $self->citation;

    for my $citation (@citations){
        $citation->base_url($self->base_url);
    }

    return 1;
}

1;
