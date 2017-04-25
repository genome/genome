package Genome::Model::Build::Command::MoveAllocations;

use strict;
use warnings;

use Genome;
use Set::Scalar;

class Genome::Model::Build::Command::MoveAllocations {
    is => 'Command::V2',
    has => [
        builds => {
            is                  => 'Genome::Model::Build',
            is_many             => 1,
            shell_args_position => 1,
            doc => 'Build(s) to use. Resolved from command line via text string.',
            require_user_verify => 0,
        },
        disk_group => {
            is                  => 'Genome::Disk::Group',
            is_optional         => 1,
            doc => 'A specific disk group to move all data to (defaults to AnP configured location)',
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

sub execute {
    my $self = shift;

    my @builds = $self->builds;
    my $total_builds = scalar @builds;
    my $current_build_num = 0;


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
        my %result_classes = map { $_->class => 1 } grep { $_->disk_allocation } @disk_usage_results;
        my $class_set_by_results = Set::Scalar->new(sort keys %result_classes);

        my $diff_associated_allocations = $class_set_by_associated_allocations->difference($class_set_by_results);
        if (my @m = $diff_associated_allocations->members) {
            $self->status_message('Skipping the following classes not defined as disk_usage_results: '. "\n". join("\n", grep { $_ ne $build->class } @m));
        }

        my $diff_results = $class_set_by_results->difference($class_set_by_associated_allocations);
        if ($diff_results->members) {
            $self->fatal_message('Found additional classes: '. "\n". join("\n", $diff_results->members));
        }

        my @allocations_to_move = grep { $result_classes{$_->owner_class_name} } @associated_allocations;
        push @allocations_to_move, $build->disk_allocation;

        my $disk_group = $self->disk_group // $self->_resolve_disk_group_for_build($build);

        local $ENV{UR_NO_REQUIRE_USER_VERIFY} = ($ENV{UR_NO_REQUIRE_USER_VERIFY} // 1);
        my $move_cmd = Genome::Disk::Command::Allocation::Move->create(
            allocations => \@allocations_to_move,
            target_group => $disk_group,
        );
        unless ($move_cmd) {
            $self->fatal_message('Failed to create move command for build %s.', $build->__display_name__);
        }
        $move_cmd->execute or
            $self->fatal_message('Failed to move allocations for build %s.', $build->__display_name__);
    }
    return 1;
}

sub _resolve_disk_group_for_build {
    my $self = shift;
    my $build = shift;

    my $anp = $build->model->analysis_project;
    unless ($anp) {
        $self->fatal_message('No analysis project defined for build %s and no disk_group specified!', $build->__display_name__);
    }

    my $config_dir = $anp->environment_config_dir;
    local $ENV{XGENOME_CONFIG_PROJECT_DIR} = $config_dir;

    my $dg = Genome::Config::get('disk_group_models');

    return Genome::Disk::Group->get(disk_group_name => $dg);
}


1;
