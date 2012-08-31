package Genome::Model::Command::SetDoNotArchive;

use strict;
use warnings;
use Genome;

class Genome::Model::Command::SetDoNotArchive {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            doc => 'models to have their latest complete build marked so they can\'t be archived',
        },
    ],
    has_optional => [
        reason => {
            is => 'Text',
            doc => 'reason provided models are to be marked',
        },
    ],
};

sub help_detail {
    return 'The last complete build for each model is marked such that their allocations can\'t be ' .
        'archived. This process marks all allocations those builds use, including instrument data, ' .
        'alignment results, variation detection results, etc';
}
sub help_brief { return 'marks last complete build of each model can unable to be archived' };
sub help_synopsis { return help_brief() . "\n" };

sub execute {
    my $self = shift;
    for my $build (map { $_->last_complete_build } $self->models) {
        next unless $build;
        for my $allocation ($build->all_allocations) {
            next unless $allocation;
            next unless $allocation->archivable;
            $allocation->archivable(0, $self->reason);
        }
    }
    $self->status_message("Last complete build of provided models marked so they can't be archived");
    return 1;
}

1;

