package Genome::Model::Build::Command::Queue;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::Build::Command::Queue {
    is => 'Genome::Command::Base',
    doc => "Queue the starting of a build.",
    has => [
        models => {
            is => 'Genome::Model',
            is_many => 1,
            doc => 'Model(s) to build. Resolved from command line via text string.',
            shell_args_position => 1,
        },
        reason => {
            is => 'Text',
            doc => 'Why the models are having a build requested',
        },
    ],

};

sub _is_hidden_in_docs { return !Genome::Sys->current_user_is_admin };

sub sub_command_sort_position { 1 }

sub help_synopsis {
    return <<EOS;
genome model build queue --reason 'new data assigned' 1234

genome model build queue --reason 'changed annotation input' somename



EOS
}

sub help_detail {
    return <<EOS;
Request that a new build for the model be scheduled by the build cron.
EOS
}

sub execute {
    my $self = shift;

    my @models = $self->models;
    map($_->build_requested(1, $self->reason), @models);

    return 1;
}

1;

