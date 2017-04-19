package Genome::Model::Build::Command::Move;

use strict;
use warnings;

use Genome;
use Set::Scalar;

class Genome::Model::Build::Command::Move {
    is => 'Genome::Disk::Command::Allocation::ActionBase',
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

sub help_brief {
    return "Move all allocations to the Analysis Project configured disk group.";
}

sub help_detail {
    return help_brief();
}

sub disk_allocation_action {
    return 'move';
}

sub _execute {
    my $self = shift;
    my @builds = $self->builds;
    my $total_builds = scalar @builds;
    my $current_build_num = 0;

    my %failed_allocations_per_build;

    $self->status_message("Moving allocations for $total_builds builds");

    BUILD: for my $build (@builds) {
        # Print progress information
        $current_build_num++;
        $self->status_message("\nWorking on build " . $build->id .
            " ($current_build_num of $total_builds)");

        my @associated_allocations = $build->associated_disk_allocations;
        my $num_allocations = @associated_allocations;
        my %owner_classes = map { $_->owner_class_name => 1 } @associated_allocations;
        my $class_set_by_associated_allocations = Set::Scalar->new(sort keys %owner_classes);


        my @disk_usage_results = $build->disk_usage_results;
        my %result_classes = map { $_->class => 1 } @disk_usage_results;
        my $class_set_by_results = Set::Scalar->new(sort keys %result_classes);

        my $diff_associated_allocations = $class_set_by_associated_allocations->difference($class_set_by_results);
        if ($diff_associated_allocations->members) {
            $self->status_message('Skipping the following classes not defined as disk_usage_results: '. "\n". join("\n",$diff_associated_allocations->members));
        }

        my $diff_results = $class_set_by_results->difference($class_set_by_associated_allocations);
        if ($diff_results->members) {
            $self->fatal_message('Found additional classes: '. "\n". join("\n", $diff_results->members));
        }
    }
    return 1;
}


1;
