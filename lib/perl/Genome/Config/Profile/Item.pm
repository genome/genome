package Genome::Config::Profile::Item;

use strict;
use warnings;

use Genome;
use File::Spec;
use File::Basename;
use Lingua::EN::Inflect;

class Genome::Config::Profile::Item {
    roles => [ 'Genome::Role::ObjectWithTimestamps', 'Genome::Role::ObjectWithCreatedBy' ],
    table_name => 'config.profile_item',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        allocation => {
            is => 'Genome::Disk::Allocation',
            reverse_as => 'owner',
            is_many => 1,
        },
        analysis_menu_item => {
            is => 'Genome::Config::AnalysisMenu::Item',
            id_by => 'analysismenu_item_id',
            is_optional => 1,
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            id_by => 'analysis_project_id',
            constraint_name => 'profile_item_analysis_project_id_fkey',
        },
        status => {
            is => 'Text',
            valid_values => [ "disabled", "active", "inactive" ],
            default_value => 'active',
        },
        model_bridges => {
            is => 'Genome::Config::AnalysisProject::ModelBridge',
            reverse_as => 'config_profile_item',
            is_many => 1,
        },
        models => {
            is => 'Genome::Model',
            via => 'model_bridges',
            to => 'model',
            is_many => 1,
        },
        tag_bridges => {
            is => 'Genome::Config::Tag::Profile::Item',
            reverse_as => 'profile_item',
            is_many => 1,
        },
        tags => {
            is => 'Genome::Config::Tag',
            via => 'tag_bridges',
            to => 'tag',
            is_many => 1,
            is_mutable => 1,
        },
        tag_names => {
            is => 'Text',
            via => 'tags',
            to => 'name',
            is_many => 1,
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

UR::Observer->register_callback(subject_class_name => __PACKAGE__, aspect => 'create', callback => \&_is_created);

sub delete {
    my $self = shift;
    eval {
        if ($self->allocation) {
            $self->allocation->delete();
        }
    };
    if(my $error = $@) {
        die($error);
    }
    return $self->SUPER::delete();
}

sub file_path {
    my $self = shift;

    if ($self->is_concrete) {
        return $self->_get_concrete_file_path();
    } else {
        return $self->analysis_menu_item->file_path;
    }
};

sub concretize {
    my $self = shift;

    return if $self->is_concrete();

    my $original_file_path = $self->file_path;
    my $skip_allocation_path_creation = 0;
    return $self->_create_allocation_for_file($original_file_path, $skip_allocation_path_creation);
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

sub has_model_for {
    my $self = shift;
    my $instrument_data = shift;

    my $model_set = Genome::Model->define_set(
        'analysis_project_bridges.profile_item_id' => $self->id,
        'instrument_data.id' => $instrument_data->id,
    );

    return $model_set->count;
}

sub _create_allocation_for_file {
    my $self = shift;
    my $file_to_store = shift;
    my $skip_allocation_path_creation = (shift // 1) && !$ENV{UR_DBI_NO_COMMIT};

    my $allocation = Genome::Disk::Allocation->create(
        owner_id            => $self->id,
        disk_group_name     => Genome::Config::get('disk_group_references'),
        allocation_path     => 'analysis_configuration/' . $self->id,
        owner_class_name    => $self->class,
        kilobytes_requested => $self->_get_size_in_kb($file_to_store),
        skip_allocation_path_creation => $skip_allocation_path_creation,
    );

    return $self->_copy_file_to_allocation($file_to_store, $allocation, $skip_allocation_path_creation);
}

sub _copy_file_to_allocation {
    my $self = shift;
    my $original_file_path = shift;
    my $allocation = shift;
    my $skip_allocation_path_creation = shift // 1;

    my $filename = File::Basename::basename($original_file_path);
    my $destination_file_path = File::Spec->catdir($allocation->absolute_path, $filename);
    if ($ENV{UR_DBI_NO_COMMIT} or not $skip_allocation_path_creation) {
        Genome::Sys->copy_file($original_file_path, $destination_file_path);
        $allocation->reallocate();
    }
    else {
        my ($file_name) = File::Basename::fileparse($original_file_path, qr/\.[^.]+/);
        my $config_file = join('.', 'model', $self->id, $file_name, 'yaml');
        Genome::Sys->shellcmd(
            cmd => [
                '/usr/bin/gsutil/gsutil',
                'cp',
                $original_file_path,
                Genome::Config::get('gcp_config_bucket') . $config_file,
            ],
            input_files => [$original_file_path],
        );
        $self->status_message('Config file queued for installation')
    }
    $allocation->archivable(0);
    return $destination_file_path;
}

sub is_concrete { return defined(shift->allocation); }

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

sub _is_created {
    my $self = shift;
    Genome::Timeline::Event::AnalysisProject->config_added(
        $self->id,
        $self->analysis_project,
    );
}

sub is_current {
    my $self = shift;

    return $self->status ne 'disabled' && $self->analysis_project->is_current;
}

1;
