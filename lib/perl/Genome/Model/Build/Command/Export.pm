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
        create_tarball => {
            is => 'Boolean',
            is_optional => 1,
            is_input => 1,
            default => 0,
            doc =>'create a .tar file containing the build directory. Useful for globus transfers that eat symlinks',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->validate_export_directory;

    my $build = $self->build;
    my $allocation = Genome::Disk::Allocation->get_allocation_for_path($build->data_directory);
    #This will not give us the full space needed but it's a start
    $self->check_available_space($allocation);

    my $export_directory = $self->export_directory;
    my $tempdir;
    if($self->create_tarball){
        #create temp directory for pre-tarball files
        $tempdir = Genome::Sys->create_temp_directory();
        unless($tempdir) {
            $self->fatal_message("Unable to create temporary directory $!");
        }
        $export_directory = File::Spec->join($tempdir, "build" . $self->build->id);
    }
   

    $allocation->copy(output_dir => $export_directory);

    my $resolve_dangling_symlinks;
    my $resolve_symlinks;
    my $find = sub {
        find({
            wanted => $resolve_symlinks,
            follow => 0,
        }, shift)
    };

    $resolve_symlinks = sub {
        my $path = $File::Find::name;
        if (-l $path) {
            my $symlink_target = abs_path($path);
            if (-e $symlink_target) {
                #Only do this if symlink doesn't point to somewhere within
                #this directory structure
                unless (index($symlink_target, $export_directory) == 0) {
                    if (index($path, $export_directory) == 0) {
                        unlink $path;
                    }
                    else {
                        Genome::Carp::confessf('Path %s escaped from export directory', $path);
                    }
                    if (-f $symlink_target) {
                        Genome::Sys->copy_file($symlink_target, $path);
                    }
                    elsif (-d $symlink_target) {
                        Genome::Sys->rsync_directory(source_directory => $symlink_target, target_directory => $path);

                        $find->($path);
                    }
                }
            }
            else {
                $resolve_dangling_symlinks->($path, $symlink_target);
            }
        }
    };

    $resolve_dangling_symlinks = sub {
        my ($path, $symlink_target) = @_;
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

    if($self->create_tarball){
        my $tarball = File::Spec->join($self->target_export_directory, "build" . $self->build->id . ".tar");   
        my $rv = Genome::Sys->shellcmd(cmd => ["tar","-cf",$tarball,"-C",$tempdir,"build" . $self->build->id]);
        unless ($rv) {
            $self->fatal_message("Could not create tar file");
        }
    };

    #TODO: Handle software result allocation that are being used by this build but aren't
    #symlinked in the build's data directory

    return 1;
}

sub validate_export_directory {
    my $self = shift;
    if (Genome::Disk::Allocation->get_allocation_for_path($self->target_export_directory)) {
        $self->fatal_message(
            "Target export directory (%s) is a path in the apipe allocation system. Please specify an external directory.",
            $self->target_export_directory
        );
    }
}

sub check_available_space {
    my $self = shift;
    my $allocation = shift;

    my $total_directory_size = Genome::Sys->directory_size_recursive($allocation->absolute_path);
    my $retirement_path = $self->target_export_directory;

    unless(-d $retirement_path) {
        $self->fatal_message('Did not find directory at <%s>.', $retirement_path);
    }

    my $disk_info = df($retirement_path, 1);
    if (!defined($disk_info) || $disk_info->{bavail} <= $total_directory_size) {
        $self->fatal_message("Not enough space available on target path (%s)", $retirement_path);
    }
}

sub export_directory {
    my $self = shift;
    return File::Spec->join($self->target_export_directory, "build" . $self->build->id);
}

1;
