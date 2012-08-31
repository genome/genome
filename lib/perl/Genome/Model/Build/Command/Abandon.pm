package Genome::Model::Build::Command::Abandon;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Abandon {
    is  => 'Genome::Model::Build::Command::Base',
    has => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc => 'Build(s) to use. Resolved from command line via text string.',
            require_user_verify => 1,
        },
    ],
};

sub sub_command_sort_position { 5 }

sub help_brief {
    return "Abandon a build and it's events";
}

sub help_detail {
    return help_brief();
}

sub execute {
    my $self = shift;

    my @builds      = $self->builds;
    my $build_count = scalar(@builds);
    for my $build (@builds) {
        $self->_total_command_count($self->_total_command_count + 1);
        my $transaction = UR::Context::Transaction->begin();
        my $successful = eval { $build->abandon };
        if ($successful and $transaction->commit) {
            $self->status_message( "Successfully abandoned build (" . $build->__display_name__ . ")." );
        }
        else {
            $self->append_error($build->__display_name__, "Failed to abandon build: $@.");
            $transaction->rollback;
        }
    }

    $self->display_command_summary_report();

    return !scalar(keys %{$self->_command_errors});
}

1;

#$HeadURL$
#$Id$
