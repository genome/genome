package Genome::Model::Build::Command::Abandon;

use strict;
use warnings;

use Genome;
use Try::Tiny qw(try catch);

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
        header_text => {
            is => 'Text',
            is_optional => 1,
            doc => 'Header for build note.',
        },
        body_text => {
            is => 'Text',
            is_optional => 1,
            doc => 'Body for build note.',
        },
    ],
    has_transient_optional => [
        show_display_command_summary_report => {
            default_value => 1,
            doc => 'whether to display the summary after abandoning builds',
        },
    ],
};

sub sub_command_sort_position { 5 }

sub help_brief {
    return "Abandon a build and its events";
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

        try {
            $build->abandon($self->header_text, $self->body_text);
            $self->successfully_abandoned_callback($build);
            $transaction->commit() or die "commit failed";
        } catch {
            $self->append_error($build->__display_name__, "Failed to abandon build: $_.");
            $transaction->rollback;
        };
    }

    $self->display_command_summary_report() if $self->show_display_command_summary_report;

    return !scalar(keys %{$self->_command_errors});
}

sub successfully_abandoned_callback {
    my $self = shift;
    my $build = shift;

    $self->status_message( "Successfully abandoned build (" . $build->__display_name__ . ")." );
}

1;

#$HeadURL$
#$Id$
