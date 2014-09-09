package Genome::Config::AnalysisProject::Command::TagConfigFile;

use strict;
use warnings;

use Genome;

class Genome::Config::AnalysisProject::Command::TagConfigFile {
    is => 'Command::V2',
    has_input => [
        tag => {
            is => 'Genome::Config::Tag',
            doc => 'The tag to apply to the configuration item(s)',
        },
        profile_items => {
            is => 'Genome::Config::Profile::Item',
            doc => 'The configuration item(s) to tag',
            shell_args_position => 1,
            is_many => 1,
        },
    ],
    doc => 'command to add a tag to individual configurations within an analysis project',
};

sub help_detail {
    return <<EOHELP
This command is used to add a tag to one or more configurations for an Analysis Project.  Once tagged, models can be filtered by the existence of these tags.

Subject mappings can also contain tags. This is used further limit what instrument data will be processed by the configuration, even if it otherwise matches.
EOHELP
}

sub help_synopsis {
    return <<EOHELP
genome config analysis-project tag-config-file --tag "my tag 1" 1234567890abcdef
EOHELP
}

sub help_brief {
    return 'tag configurations in an analysis-project'
}

sub execute {
    my $self = shift;

    my $tag = $self->tag;
    my @items = $self->profile_items;
    my @bridges;
    for my $item (@items) {
        push @bridges, Genome::Config::Tag::Profile::Item->create(
            tag => $tag,
            profile_item => $item,
        );
    }

    $self->status_message('Configurations tagged: %d', scalar(@bridges));
    return 1;
}

1;
