package Genome::Site::TGI::Measurable; 

use strict;
use warnings;

use Genome;

require Carp;
require File::Basename;
require File::Copy;

class Genome::Site::TGI::Measurable {
    table_name => 'PHENOTYPE_MEASURABLE',
    is_abstract => 1,
    subclassify_by => '_subclass_by_subject_type',
    id_by => [
        subject_id => { 
            is => 'Number',
            doc => 'the numeric ID for the specimen in both the LIMS and the analysis system', 
        },
    ],
    has => [
        subject_type => {
            is => 'Text',
            column_name =>'SUBJECT_TYPE',
        },
        # These are here, and should be overidden in the subclass
        name => { column_name => '', },
        common_name => { calculate => q| return $_[0]->name; |, },
        # disk
        disk_allocation => {
            is => 'Genome::Disk::Allocation', 
            calculate_from => [ 'class', 'id' ],
            calculate => sub{
                my ($class, $id) = @_;
                my $disk_allocation = Genome::Disk::Allocation->get(
                    owner_class_name => $class,
                    owner_id => $id,
                );
                return $disk_allocation;
            },
        },
        data_directory => { 
            is => 'Text', 
            calculate_from => [qw/ disk_allocation /],
            calculate => sub{
                my $disk_allocation = shift;
                return if not $disk_allocation;
                return $disk_allocation->absolute_path;
            },
            #is => 'Text', via => 'disk_allocation', to => 'absolute_path', },
        },
    ],
    data_source => 'Genome::DataSource::Dwrac',
};

sub __display_name__ {
    return $_[0]->name.' ('.$_[0]->id.')';
}

sub _subclass_by_subject_type {
    my ($measurable) = @_;

    #print Data::Dumper::Dumper(\@_);
    my $subject_type = $measurable->subject_type;
    if ( $subject_type eq 'organism sample' or $subject_type eq 'organism_sample' ) {
        return 'Genome::Site::TGI::Sample';
    }
    elsif ( $subject_type eq 'organism taxon' or $subject_type eq 'organism_taxon' ) {
        return 'Genome::Site::TGI::Taxon';
    }
    elsif ( $subject_type eq 'organism individual' or $subject_type eq 'organism_individual' ) {
        return 'Genome::Site::TGI::Individual';
    }
    elsif ( $subject_type eq 'population group' or $subject_type eq 'population_group' ) {
        return 'Genome::Site::TGI::Individual';
    }
    else {
        Carp::confess("Unknown subject type ($subject_type), can't determine approporate subclass");
    }
}

sub add_file {
    my ($self, $file) = @_;

    $self->status_message('Add file to '. $self->__display_name__);

    Carp::confess('No file to add') if not $file;
    my $size = -s $file;
    Carp::confess("File ($file) to add does not have any size") if not $size;
    my $base_name = File::Basename::basename($file);
    Carp::confess("Could not get basename for file ($file)") if not $base_name;
    
    my $disk_allocation = $self->disk_allocation;
    if ( not $disk_allocation ) {
        # Create
        $disk_allocation = Genome::Disk::Allocation->allocate(
            disk_group_name => 'info_genome_models',
            allocation_path => '/model_data/'.$self->id,
            kilobytes_requested => $size,
            owner_class_name => $self->class,
            owner_id => $self->id
        );
        if ( not $disk_allocation ) {
            Carp::confess('Failed to create disk allocation to add file');
        }
    }
    else { 
        # Make sure we don't overwrite
        if ( grep { $base_name eq $_ } map { File::Basename::basename($_) } glob($disk_allocation->absolute_path.'/*') ) {
            Carp::confess("File ($base_name) to add already exists in path (".$disk_allocation->absolute_path.")");
        }
        # Reallocate w/ move to accomodate the file
        my $realloc = eval{
            $disk_allocation->reallocate(
                kilobytes_requested => $disk_allocation->kilobytes_requested + $size,
                allow_reallocate_with_move => 1,
            );
        };
        if ( not $realloc ) {
            Carp::confess("Cannot reallocate (".$disk_allocation->id.") to accomadate the file ($file)");
        }
    }

    my $absolute_path = $disk_allocation->absolute_path;
    if ( not -d $absolute_path ){
        Carp::confess('Absolute path does not exist for disk allocation: '.$disk_allocation->id);
    }
    my $to = $absolute_path.'/'.$base_name;
    
    $self->status_message("Copy $file to $to");
    my $copy = File::Copy::copy($file, $to);
    if ( not $copy ) {
        Carp::confess('Copy of file failed');
    }

    my $new_size = -s $to;
    if ( $new_size != $size ) {
        Carp::confess("Copy of file ($file) succeeded, but file ($to) has different size.");
    }

    $self->status_message('Add file...OK');

    return 1;
}

sub get_files {
    my $self = shift;

    my $disk_allocation = $self->disk_allocation;
    return if not $disk_allocation;

    return glob($disk_allocation->absolute_path.'/*');
}

1;

