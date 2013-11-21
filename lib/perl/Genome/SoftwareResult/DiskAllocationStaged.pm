package Genome::SoftwareResult::DiskAllocationStaged;

use Genome;
use strict;
use warnings;
use File::Path 'remove_tree';
use Carp;
use YAML;

# This class takes care of getting a disk allocation and using a subdirectory
# to stage the result.  When the result has been completely generated the
# data are moved from the staging subdirectory to the base directory of the
# disk allocation.
#
# things you MUST override:
#  _generate_result
#
# things you might want to override:
# _allocation_subdirectory
# _allocation_disk_group_name
# _staging_kilobytes_requested

class Genome::SoftwareResult::DiskAllocationStaged {
    is => 'Genome::SoftwareResult',
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    my $disk_allocation = $self->_acquire_disk_allocation(
            $self->_allocation_subdirectory,
            $self->_allocation_disk_group_name,
            $self->_staging_kilobytes_requested);
    my $allocation_path = $disk_allocation->absolute_path;
    $self->output_dir($allocation_path);

    my $staging_directory;
    eval {
        $staging_directory = $self->_setup_staging_directory($allocation_path);
        $self->_generate_result($staging_directory);
    };
    if($@) {
        my $error_message = $@;
        my $self_error_message = $self->error_message("Result generation failed...\n$error_message");
        $self->delete();
        Carp::croak($self_error_message);
    } else {
        $self->_destage_result($staging_directory, $allocation_path);
    }

    $disk_allocation->reallocate();
    return $self;
}

# override to create the result
sub _generate_results {
    my ($self, $staging_directory) = @_;
    Carp::croak($self->error_message('_generate_results not implemented'));
}

# override to customize how disk allocation is acquired
sub _allocation_subdirectory {
    my ($self) = @_;
    # NOTE: we (apipe) currently don't have write access to any subdirectories
    # under info other than 'model_data' and 'build_merged_alignments' so...
    #                     |<------------ this must be one of those
    return join('/', 'model_data', 'software-result', $self->id);
}

# override to customize how disk allocation is acquired
sub _allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

# override to customize how disk allocation is acquired
#   How much disk will you require for staging your result?
sub _staging_kilobytes_requested {
    my ($self) = @_;
    return '5000000';
}

sub _acquire_disk_allocation {
    my ($self, $subdirectory, $disk_group_name, $kilobytes_requested) = @_;

    my %allocation_params = (
            disk_group_name => $disk_group_name,
            allocation_path => $subdirectory,
            kilobytes_requested => $kilobytes_requested,
            owner_class_name => $self->class,
            owner_id => $self->id);

    my $allocation = Genome::Disk::Allocation->allocate(
        %allocation_params);
    unless ($allocation) {
        Carp::croak($self->error_message(
                "Failed to get disk allocation with params:\n" .
                YAML::Dump(%allocation_params)));
    }

    my $allocation_path = $allocation->absolute_path;
    unless (-d $allocation_path) {
        $allocation->delete;
        Carp::croak($self->error_message(
                "Allocation path $allocation_path doesn't exist!"));
    }

    return $allocation;
}

sub _setup_staging_directory {
    my ($self, $allocation_path) = @_;
    my $staging_directory = Genome::Sys->create_directory(
            join('/', $allocation_path, 'staging_area'));
    return $staging_directory;
}

sub _delete_staging_directory {
    my ($self, $staging_directory) = @_;
    remove_tree($staging_directory);
    if(-e $staging_directory) {
        self->error_message("Couldn't remove staging_directory: $staging_directory");
    }
}

sub _destage_result {
    my ($self, $staging_directory, $allocation_path) = @_;

    my @names = Genome::Sys->list_directory($staging_directory);
    for my $name (@names) {
        my $source = join('/', $staging_directory, $name);
        my $destination = join('/', $allocation_path, $name);
        if(not -e $destination) {
            rename($source, $destination);
        } else {
            $self->error_message("Couldn't move $source to $destination since desination already exists.");
        }
    }
    $self->_delete_staging_directory($staging_directory);
}

1;
