package Genome::Model::Build::Command::Unarchive;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Unarchive {
    is  => 'Genome::Model::Build::Command::Base',
    has => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc => 'Build(s) to use. Resolved from command line via text string.',
            require_user_verify => 0,
        },
    ],
};

sub sub_command_sort_position { 7 }

sub help_brief {
    return "Unarchive all allocations related to a build";
}

sub help_detail {
    return help_brief();
}

sub _limit_results_for_builds {
    my ( $class, @builds ) = @_;
    return @builds;
}

sub execute {
    my $self = shift;

    my @builds      = $self->builds;
    my $build_count = scalar(@builds);
    for my $build (@builds) {

        my $build_id = $build->id;
        $self->_total_command_count($self->_total_command_count + 1);

        my ($alloc_count, $unarchived_count);

        my $transaction = UR::Context::Transaction->begin();
        my $successful = eval {
            my @allocations = $build->all_allocations();

            $alloc_count = scalar @allocations;
            print "Found $alloc_count allocations to unarchive...\n";
            ALLOC:
            for my $alloc (@allocations) {
                next ALLOC if ! $alloc->is_archived();

                print "unarchiving ($build_id - " . $alloc->id . "): " . $alloc->allocation_path . "\n";
                if ($alloc->unarchive()) {
                    $unarchived_count++;
                } else {
                    $self->append_error($build->__display_name__, "ERROR: Failed to unarchive alloc id : $@.");
                }
            }
            return 1;
        };

        if ($successful and $transaction->commit) {
            $self->status_message("Unarchived $unarchived_count of $alloc_count allocations for build (" 
            . $build->__display_name__ 
            . ").");
        } else {
            $self->append_error($build->__display_name__, "ERROR: success = $successful, trying to rollback");
            $transaction->rollback;
        }
    }

    $self->display_command_summary_report();

    return !scalar(keys %{$self->_command_errors});
}

1;

#$HeadURL$
#$Id$
