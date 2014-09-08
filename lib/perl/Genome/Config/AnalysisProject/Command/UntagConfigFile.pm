package Genome::Config::AnalysisProject::Command::UntagConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::UntagConfigFile {
    is => 'Command::V2',
    has_input => [
        tag => {
            is => 'Genome::Config::Tag',
            doc => 'The tag to remove from the configuration item(s)',
        },
        profile_items => {
            is => 'Genome::Config::Profile::Item',
            doc => 'The configuration item(s) to untag',
            shell_args_position => 1,
            is_many => 1,
        },
    ],
    doc => 'command to remove a tag from individual configurations within an analysis project',
};

sub help_detail {
    return <<EOHELP
This command is used to remove a tag from one or more configurations for an Analysis Project.
EOHELP
}

sub help_synopsis {
    return <<EOHELP
genome config analysis-project untag-config-file --tag "my tag 1" 1234567890abcdef
EOHELP
}

sub help_brief {
    return 'untag configurations in an analysis-project'
}

sub execute {
    my $self = shift;

    my $tag = $self->tag;
    my @items = $self->profile_items;
    my $count = 0;
    for my $item (@items) {
        my $bridge = Genome::Config::Tag::Profile::Item->get(
            tag => $tag,
            profile_item => $item,
        );
        if($bridge) {
            $bridge->delete;
            $count++;
        } else {
            $self->warning_message(
                'Tag "%s" not found on profile item %s.',
                $tag->name,
                $item->id,
            );
        }
    }

    $self->status_message('Configurations untagged: %d', $count);
    return 1;
}

1;
