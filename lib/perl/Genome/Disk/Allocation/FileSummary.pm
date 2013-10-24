package Genome::Disk::Allocation::FileSummary;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Disk::Allocation::FileSummary {
    is => ['Genome::Utility::ObjectWithTimestamps', 'Genome::Utility::ObjectWithCreatedBy'],
    table_name => 'disk.file_summary',
    id_generator => '-uuid',
    data_source => 'Genome::DataSource::GMSchema',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        }
    ],
    has => [
        allocation_id => {
            is => 'Text',
            len => 64,
        },
        allocation => {
            is => 'Genome::Disk::Allocation',
            id_by => 'allocation_id',
        },
        file => {
            is => 'Text',
            doc => 'relative path to file from the allocation root',
        },
        digest => {
            is => 'Text',
        },
        size_in_bytes => {
            is => 'Number',
        },
        is_symlink => {
            is => 'Boolean',
        },
        destination => {
            is => 'Text',
            is_optional => 1,
        }

    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return $self->initialize();
}

sub create_or_update {
    my $class = shift;
    my $self = $class->get(@_);
    unless ($self) {
       $self = $class->SUPER::create(@_);
    }
    return $self->initialize()
}

sub initialize {
    my $self = shift;

    eval {
        my $file_path = File::Spec->join($self->allocation->absolute_path, $self->file);

        if (-f $file_path) {
            $self->digest($self->get_file_digest($file_path));
            $self->size_in_bytes($self->get_file_size($file_path));
        }

        if (-l $file_path) {
            $self->is_symlink(1);
            $self->destination(readlink($file_path));
        } elsif (not -e $file_path) {
            die("Sorry, but $file_path does not exist!");
        } elsif (-d $file_path) {
            die("Sorry, but $file_path is a directory!");
        } else {
            $self->is_symlink(0);
        }
    };

    if(my $error = $@) {
        die($error);
    }

    return $self;
}

sub get_file_digest {
    my $self = shift;
    my $file_path = shift;

    my $digest = Genome::Sys->md5sum($file_path);
    die("Unable to calculate md5 for $file_path!") unless $digest;

    return $digest;
}

sub get_file_size {
    my $self = shift;
    my $file_path = shift;

    my $size = -s $file_path;
    if (defined($size)) {
        return $size;
    } else {
        die("Unable to get size for $file_path!");
    }
}

1;
