package Genome::Model::Build::Command::Export;

use strict;
use warnings;
use Genome;
use File::Find;
use Filesys::Df;
use Cwd qw(abs_path);

class Genome::Model::Build::Command::Export {
    is => 'Command::V2',
    has => [
        build => {
            is => 'Genome::Model::Build',
        },
        target_export_directory => {
            is => 'Text',
        },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $allocation = Genome::Disk::Allocation->get_allocation_for_path($build->data_directory);
    #This will not give us the full space needed but it's a start
    $self->check_available_space($allocation);

    my $export_directory = $self->export_directory;
    $allocation->copy(output_dir => $export_directory);

    my $resolve_dangling_symlinks;
    my $resolve_symlinks;
    my $find = sub {
        find({
            wanted => $resolve_symlinks,
            follow => 1,
            dangling_symlinks => $resolve_dangling_symlinks,
        }, shift)
    };

    $resolve_symlinks = sub {
        my $path = $File::Find::name;
        if (-l $path) {
            my $symlink_target = $File::Find::fullname;
            my $allocation = Genome::Disk::Allocation->get_allocation_for_path($symlink_target);
            if (defined($allocation) && $allocation->is_active) {
                unlink $path;
                $allocation->copy(output_dir => $path);
            }
            elsif (-e $symlink_target) {
                #Only do this if symlink doesn't point to somewhere within
                #this directory structure
                unless (index($symlink_target, $export_directory) == 0) {
                    unlink $path;
                    if (-f $symlink_target) {
                        Genome::Sys->copy_file($symlink_target, $path);
                    }
                    elsif (-d $symlink_target) {
                        Genome::Sys->rsync_directory(source_directory => $symlink_target, target_directory => $path);
                    }
                }
            }
            $find->($path);
        }
    };

    $resolve_dangling_symlinks = sub {
        my $path = File::Spec->join($_[1], $_[0]);
        my $symlink_target = abs_path($path);
        my $allocation = Genome::Disk::Allocation->get_allocation_for_path($symlink_target);
        unlink $path;
        if (defined($allocation)) {
            if ($allocation->is_archived || $allocation->is_active) {
                $allocation->copy(output_dir => $path);
                $find->($path);
            }
            elsif ($allocation->is_purged) {
                my $text = sprintf("Allocation (%s) has been purged and its data cannot be recovered.\n", $allocation->id);
                my $owner = $allocation->owner;
                if ($owner->isa('Genome::SoftwareResult')) {
                    $text .= sprintf("This allocation belongs to software result (%s) of class (%s) with test name (%s).\n", $owner->id, $owner->class, $owner->test_name);
                }
                Genome::Sys->write_file($path, $text);
            }
        }
        else {
            my $text = sprintf("Symlink target (%s) does not exist.\n", $symlink_target);
            Genome::Sys->write_file($path, $text);
        }
    };

    $find->($export_directory);

    #TODO: Handle software result allocation that are being used by this build but aren't
    #symlinked in the build's data directory

    return 1;
}

sub check_available_space {
    my $self = shift;
    my $allocation = shift;

    my $total_directory_size = Genome::Sys->directory_size_recursive($allocation->absolute_path);
    my $retirement_path = $self->target_export_directory;
    my $disk_info = df($retirement_path, 1);
    if (!defined($disk_info) || $disk_info->{bavail} <= $total_directory_size) {
        $self->fatal_message("Not enough space available on target path (%s)", $retirement_path);
    }
}

sub export_directory {
    my $self = shift;
    return File::Spec->join($self->target_export_directory, $self->build->id);
}

1;
