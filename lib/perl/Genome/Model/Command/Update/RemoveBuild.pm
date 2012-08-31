package Genome::Model::Command::Update::RemoveBuild;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Update::RemoveBuild {
    is => 'Genome::Command::Base',
    has => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc => 'Build(s) to use. Resolved from command line via text string.',
            require_user_verify => 1,
        },
    ],
    has_optional => [
        keep_build_directory => {
            is => 'Boolean',
            default_value => 0,
            doc => 'A boolean flag to allow the retention of the model directory after the model is purged from the database.(default_value=0)',
        },
    ],
};

sub sub_command_sort_position { 7 }

sub _is_hidden_in_docs { return !Genome::Sys->current_user_is_admin };

sub help_brief {
    "Remove a build.";
}

sub help_detail {
    "This command will remove a build from the system.  The rest of the model remains the same, as does independent data like alignments.";
}

# TODO This needs to be cleaned up
sub execute {
    my $self = shift;

    unless (Genome::Sys->current_user_is_admin) {
        print "You must be an admin to use this command.\n";
        return;
    }

    my @builds = $self->builds;
    my $build_count = scalar(@builds);
    for my $build (@builds) {
        $self->_total_command_count($self->_total_command_count + 1);
        my $transaction = UR::Context::Transaction->begin();
        my $display_name = $build->__display_name__;
        my $successful = eval {
            my @__errors__ = $build->__errors__;
            if (@__errors__) {
                die "build or instrument data has __errors__, cannot remove: " . join('; ', @__errors__);
            }
            $build->delete;
        };
        if ($successful and $transaction->commit) {
            $self->status_message("Successfully removed build (" . $display_name . ").");
        }
        else {
            $self->append_error($display_name, "Failed to remove build: $@.");
            $transaction->rollback();
        }
    }

    $self->display_command_summary_report();

    return !scalar(keys %{$self->_command_errors});
}

1;
