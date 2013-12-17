package Genome::Config::Profile::Item;

use strict;
use warnings;

use Genome;
use File::Spec;
use File::Basename;
use Lingua::EN::Inflect;

class Genome::Config::Profile::Item {
    is => ['Genome::Utility::ObjectWithTimestamps','Genome::Utility::ObjectWithCreatedBy'],
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'config.profile_item',
    id_generator => '-uuid',
    has => [
        id => {
            is => 'Text',
            len => 64,
        },
        allocation => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            reverse_as => 'owner'
        },
        analysis_menu_item => {
            is => 'Genome::Config::AnalysisMenu::Item',
            id_by => 'analysismenu_item_id',
            is_optional => 1,
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
        },
    ],
};

sub file_path {
    my $self = shift;

    if ($self->_is_concrete) {
        return $self->_get_concrete_file_path();
    } else {
        return $self->analysis_menu_item->file_path;
    }
};

sub concretize {
    my $self = shift;

    return if $self->_is_concrete();

    my $original_file_path = $self->file_path;
    return $self->_create_allocation_for_file($original_file_path);
}

sub create_from_file_path {
    my $class = shift;
    my %params = @_;
    my $file = delete $params{file_path};

    die('Must supply the path to a file as "file_path"!') unless $file;

    my $profile_item = $class->create(%params);
    $profile_item->_create_allocation_for_file($file);
    return $profile_item;
}

sub _create_allocation_for_file {
    my $self = shift;
    my $file_to_store = shift;

    my $allocation = Genome::Disk::Allocation->create(
        owner_id            => $self->id,
        disk_group_name     => 'info_apipe_ref',
        allocation_path     => 'analysis_configuration/' . $self->id,
        owner_class_name    => $self->class,
        kilobytes_requested => $self->_get_size_in_kb($file_to_store),
    );

    return $self->_copy_file_to_allocation($file_to_store, $allocation);
}

sub _copy_file_to_allocation {
    my $self = shift;
    my $original_file_path = shift;
    my $allocation = shift;

    my $filename = File::Basename::basename($original_file_path);
    my $destination_file_path = File::Spec->catdir($allocation->absolute_path, $filename);
    Genome::Sys->copy_file($original_file_path, $destination_file_path);
    $allocation->reallocate();
    return $destination_file_path;
}

sub _is_concrete { return defined(shift->allocation); }

sub _get_concrete_file_path {
    my $self = shift;

    my @files = glob($self->allocation->absolute_path . '/*');
    if (scalar(@files) != 1) {
        die(sprintf('%s found when one expected!',
            Lingua::EN::Inflect::NO('file', scalar(@files))));
    }

    return $files[0];
}

sub _get_size_in_kb { return ((-s $_[1])/1024) + 1; }

1;
