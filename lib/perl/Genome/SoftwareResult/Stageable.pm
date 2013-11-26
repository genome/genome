package Genome::SoftwareResult::Stageable;

use warnings;
use strict;
use Genome;
use Sys::Hostname;
use File::Path;

class Genome::SoftwareResult::Stageable {
    is => 'Genome::SoftwareResult',
    is_abstract => 1,
    has_transient => [
        temp_staging_directory => {
            is => 'Text',
            doc => 'Directory to use for staging the generated data before putting on allocated disk.',
            is_optional => 1
        }

    ]
};

sub resolve_allocation_kilobytes_requested {
    return $_[0]->_staging_disk_usage;
}

sub resolve_allocation_subdirectory {
    die "Must define resolve_allocation_subdirectory in your subclass of Genome::SoftwareResult::Stageable";
}

sub resolve_allocation_disk_group_name {
    die "Must define resolve_allocation_disk_group_name in your subclass of Genome::SoftwareResult::Stageable";
}

sub _working_dir_prefix {
    "software-result";
}

sub _prepare_staging_directory {
    my $self = shift;

    return $self->temp_staging_directory if ($self->temp_staging_directory);

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $basedir = sprintf("%s-%s-%s-%s-%s", $self->_working_dir_prefix, $hostname, $user, $$, $self->id);
    my $tempdir = Genome::Sys->create_temp_directory($basedir);
    unless(-d $tempdir) {
        die "failed to create a temp staging directory for completed files";
    }
    $self->temp_staging_directory($tempdir);


    return $tempdir;
}

sub _staging_disk_usage {

    my $self = shift;
    my $usage;
    unless (defined($usage = Genome::Sys->disk_usage_for_path($self->temp_staging_directory))) {
        $self->error_message("Failed to get disk usage for staging: " . Genome::Sys->error_message);
        die $self->error_message;
    }

    return $usage;
}

sub _needs_symlinks_followed_when_syncing {
    return 0;
}


sub _promote_data {
    my $self = shift;

    #my $container_dir = File::Basename::dirname($self->output_dir);
    my $staging_dir = $self->temp_staging_directory;
    my $output_dir  = $self->output_dir;

    unless ($output_dir) {
        die $self->error_message("output_dir not set");
    }

    unless ($staging_dir) {
        die $self->error_message("staging_dir not set");
    }

    $self->status_message("Now de-staging data from $staging_dir into $output_dir"); 

    my $cp_params = 'r';
    $cp_params .= 'L' if ($self->_needs_symlinks_followed_when_syncing);
    my $copy_cmd = sprintf("cp -$cp_params %s/* %s/", $staging_dir, $output_dir);
    $self->status_message("Running cp: $copy_cmd");
    my $copy_exit_code = system($copy_cmd);

    if ($copy_exit_code != 0) {

        my $rsync_params = "-avz";
        $rsync_params .= "L" if ($self->_needs_symlinks_followed_when_syncing);

        my $rsync_cmd = sprintf("rsync %s %s/* %s/", $rsync_params, $staging_dir, $output_dir);

        my $rsync_exit_code = system($rsync_cmd);
        $self->status_message("Running Rsync: $rsync_cmd");

        unless ($rsync_exit_code == 0) {
            $self->error_message("Did not get a valid return from rsync, exit code was $rsync_exit_code for call $rsync_cmd.  Cleaning up and bailing out");
            rmtree($output_dir);
            my @a = $self->disk_allocations;
            for my $a (@a) {
                $a->delete;
            }
            die $self->error_message;
        }
    }

    chmod 02775, $output_dir;
    for my $subdir (grep { -d $_  } glob("$output_dir/*")) {
        chmod 02775, $subdir;
    }

    # Make everything in here read-only 
    for my $file (grep { -f $_  } glob("$output_dir/*")) {
        chmod 0444, $file;
    }

    $self->status_message("Files in $output_dir: \n" . join "\n", glob($output_dir . "/*"));

    return $output_dir;
}

sub _reallocate_disk_allocation {
    my $self = shift;
    my $allocation = $self->disk_allocations;
    unless ($allocation) {
        $self->status_message("No allocations to resize/reallocate.");
        return 1;
    }
    $self->status_message('Resizing the disk allocation...');
    my $rv = eval { $allocation->reallocate };
    my $error = $@;
    if (!$rv) {
        my $warning_message = 'Failed to reallocate disk allocation (' . $allocation->__display_name__ . ').';
        $warning_message   .= " Error: '$error'." if $error;
        $self->warning_message($warning_message);
    }
    return $rv;
}

1;
