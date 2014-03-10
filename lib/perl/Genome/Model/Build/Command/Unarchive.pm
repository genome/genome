package Genome::Model::Build::Command::Unarchive;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Unarchive {
    is => 'Genome::Disk::Command::Allocation::UnarchiveBase',
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
    return "Unarchive all allocations related to a build";
}

sub help_detail {
    return help_brief();
}

sub _execute {
    my $self = shift;
    my @builds = $self->builds;
    my $total_builds = scalar @builds;
    my $current_build_num = 0;

    my %failed_allocations_per_build;

    # Create a temp directory on network accessible disk for logs.
    # If no errors occur, this is cleaned up. Otherwise the user is 
    # told to clean it up after reviewing any failures.
    my $dir = File::Temp::tempdir(
        'build_unarchive_XXXXXX', 
        DIR => "/gsc/var/tmp/build_unarchive", 
        CLEANUP => 0
    );
    chmod 0755, $dir;

    $self->status_message("Unarchiving data for $total_builds builds, logs being written to $dir");

    BUILD: for my $build (@builds) {
        # Print progress information
        $current_build_num++;
        $self->status_message("\nWorking on build " . $build->id . 
            " ($current_build_num of $total_builds)");

        my @allocations = grep { $_->is_archived } $build->all_allocations;
        my $num_allocations = @allocations;
        my (%job_statuses, %job_to_allocation_mapping);
        if ( @allocations ) {
            $self->status_message("Found $num_allocations archived allocations related to this build, " .
                "starting unarchive process now.");
            # Schedule unarchive of all allocations related to build via LSF
            my %jobs_to_allocations = $self->_bsub_unarchives($dir, @allocations);
            %job_to_allocation_mapping = %jobs_to_allocations;
            my @job_ids = keys %job_to_allocation_mapping;
            $self->debug_message("Unarchives scheduled, waiting for completion");
            # Wait for all jobs to finish and gather status
            my %statuses = Genome::Sys->wait_for_lsf_jobs(@job_ids);
            %job_statuses = (%job_statuses, %statuses);

            # since unarchives are happening in another process we need to reload them
            my $ids = [map { $_->id } @allocations];
            UR::Context->current->reload('Genome::Disk::Allocation', id => $ids);
        }

        # Old builds may have allocations symlinked to the data dir, and they may are archived, or the link is broken
        my @symlinked_allocations_that_need_unarchiving = grep { $_->is_archived } $build->symlinked_allocations;
        if ( @symlinked_allocations_that_need_unarchiving ) {
            $num_allocations += @symlinked_allocations_that_need_unarchiving;
            $self->status_message("Found ".@symlinked_allocations_that_need_unarchiving." archived symlinked allocations. Unarchiving...");
            my %jobs_to_allocations = $self->_bsub_unarchives($dir, @symlinked_allocations_that_need_unarchiving);
            %job_to_allocation_mapping = %jobs_to_allocations;
            my @symlinked_allocation_job_ids = keys %jobs_to_allocations;
            $self->debug_message("Unarchives for symlinked allocations scheduled, waiting for completion");
            my %symlinked_allocation_job_statuses = Genome::Sys->wait_for_lsf_jobs(@symlinked_allocation_job_ids);
            %job_statuses = (%job_statuses, %symlinked_allocation_job_statuses);

            # since unarchives are happening in another process we need to reload them
            my $ids = [map { $_->id } @symlinked_allocations_that_need_unarchiving];
            UR::Context->current->reload('Genome::Disk::Allocation', id => $ids);
        }


        # Relink broken symlinked allocations
        $build->relink_symlinked_allocations();
        my @input_builds = $build->input_builds;
        $self->status_message(sprintf("Found %s input-builds related to this build (%s).",
                scalar(@input_builds), $build->id));
        for my $input_build (@input_builds) {
            $self->status_message(sprintf("Relinking symlinked allocations on build %s",
                    $input_build->id));
            $input_build->relink_symlinked_allocations();
        }

        # Set the data directory to the absolute path of the main alloc
        my $main_allocation = $build->disk_allocation;
        UR::Context->reload($main_allocation); # this may have been loaded and unarchived earlier
        if ( $build->data_directory ne $main_allocation->absolute_path ) {
            $build->data_directory($main_allocation->absolute_path);
        }

        if ( $num_allocations == 0 ) {
            $self->status_message('No unarchived allocations found!');
            next;
        }

        # Find failures and report them (now and once all builds have completed)
        $self->debug_message("Unarchived finished, combing over return statuses and gathering up any failures");
        my %failed_allocations = $self->_get_failed_unarchives(\%job_statuses, \%job_to_allocation_mapping, $dir);
        my $fail_count = %failed_allocations ? keys %failed_allocations : 0;
        $failed_allocations_per_build{$build->id} = \%failed_allocations if %failed_allocations;

        # Print failure summary for this batch
        my $msg = "Finished build " . $build->id;
        if ($fail_count > 0) {
            $msg .= ", $fail_count of $num_allocations failed. All failures will be printed in a " .
            "summary report before this command exits.";
        }
        else {
            $msg .= ", all $num_allocations unarchives finished successfully.";
        }
        $self->status_message($msg);
    }

    # Print error summary (if needed) and exit
    $self->status_message("\nUnarchives for all builds have completed.");
    if (%failed_allocations_per_build) {
        $self->status_message("Errors occurred during unarchiving. An error summary follows. " .
            "Once all errors have been addressed, please remove the log directory $dir");
        $self->_print_error_summary(%failed_allocations_per_build);
        die "Encountered errors during unarchiving!";
    }
    else {
        $self->status_message("No errors found while unarchiving, removing log directory $dir!");
        $self->_remove_log_dir($dir);
    }

    return 1;
}

sub _bsub_unarchives {
    my $self = shift;
    my $log_file_dir = shift;
    my @allocations = @_;

    my $lab = $self->lab;
    my $requestor = $self->requestor->id;

    my %job_to_allocation_mapping;
    for my $allocation (@allocations) {
        my $allocation_id = $allocation->id;
        my @cmd = ('genome', 'disk', 'allocation', 'unarchive', $allocation_id,
            '--lab', $lab, '--requestor', $requestor,
        );
        my $job_id = Genome::Sys->bsub(
            queue => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
            cmd => \@cmd,
            log_file => "$log_file_dir/$allocation_id.out",
            err_file => "$log_file_dir/$allocation_id.err",
            job_group => '/apipe/build-unarchive',
        );
        $job_to_allocation_mapping{$job_id} = $allocation_id;
        $self->debug_message("Scheduled unarchive of allocation $allocation_id  " .
            "via LSF job $job_id");
    }
    return %job_to_allocation_mapping;
}

sub _get_failed_unarchives {
    my $self = shift;
    my %job_statuses = %{shift @_};
    my %job_to_allocation_mapping = %{shift @_};
    my $dir = shift;

    my %failed_allocations;
    for my $job_id (sort keys %job_statuses) {
        my $lsf_status = $job_statuses{$job_id};
        next unless $lsf_status eq 'EXIT';

        my $allocation_id = $job_to_allocation_mapping{$job_id};
        $self->warning_message("Allocation $allocation_id failed to unarchive via LSF job $job_id");

        $failed_allocations{$allocation_id} = [$job_id, "$dir/$allocation_id"];
    }
    return %failed_allocations;
}

sub _remove_log_dir {
    my ($self, $dir) = @_;
    my $rv = eval { Genome::Sys->remove_directory_tree($dir) };
    my $error = $@;
    if ($error) {
        $self->warning_message("Could not remove log directory $dir, please fix and remove manually.");
    }
    else {
        $self->status_message("Log directory $dir removed!");
    }
    return 1;
}

sub _print_error_summary {
    my $self = shift;
    my %failed_allocations_per_build = @_;

    $self->status_message("Error summary:");
    my @lines;
    push @lines, join("\t", 'allocation_id', 'build_id', 'lsf_job_id', 'log_file');

    for my $build (sort keys %failed_allocations_per_build) {
        my %failed_unarchives = %{$failed_allocations_per_build{$build}};
        for my $allocation (sort keys %failed_unarchives) {
            my ($job_id, $log_file) = @{$failed_unarchives{$allocation}};
            push @lines, join("\t", $allocation, $build, $job_id, $log_file);
        }
    }

    $self->status_message(join("\n", @lines));
    return 1;
}

1;
