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
    my $filename = File::Basename::basename($original_file_path);

    my $allocation = Genome::Disk::Allocation->create(
        owner_id            => $self->id,
        disk_group_name     => 'info_apipe_ref',
        allocation_path     => 'analysis_configuration/' . $self->id,
        owner_class_name    => $self->class,
        kilobytes_requested => $self->_get_allocaton_size(),
    );

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

sub _get_allocaton_size { return ((-s shift->file_path)/1024) + 1; }

1;
