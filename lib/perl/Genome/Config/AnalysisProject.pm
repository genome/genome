package Genome::Config::AnalysisProject;

use strict;
use warnings;

use Genome;

use List::MoreUtils qw(any);
use File::Spec;

class Genome::Config::AnalysisProject {
    roles => [qw(
        Genome::Role::ObjectWithTimestamps
        Genome::Role::ObjectWithCreatedBy
        Genome::Role::SoftwareResultSponsor
        Genome::Role::Searchable
        Genome::Role::ObjectWithAllocations
    )],
    table_name => 'config.analysis_project',
    id_by => [
        id => { is => 'Text', len => 64 },
    ],
    has => [
        name => { is => 'Text' },
        status => {
            is => 'Text',
            len => 255,
            default_value => 'Pending',
            valid_values => [ 'Pending', 'In Progress', 'Completed', 'Archived', 'Hold', 'Template', 'Deprecated' ],
        },
        is_cle => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Is this an analysis project for the CLIA Licensed Environment?',
        },
        run_as => {
            is => 'Text',
            len => 64,
            is_optional => 1,
            doc => 'The user account that will be used to run these models',
        },
        subject_mappings => {
            is => 'Genome::Config::AnalysisProject::SubjectMapping',
            reverse_as => 'analysis_project',
            is_many => 1,
        },
        analysis_project_bridges => {
            is => 'Genome::Config::AnalysisProject::InstrumentDataBridge',
            reverse_as => 'analysis_project',
            is_many => 1,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            via => 'analysis_project_bridges',
            is_many => 1,
        },
        samples => {
            is => 'Genome::Subject',
            via => 'instrument_data',
            to => 'sample',
            is_many => 0,
        },
        model_bridges => {
            is => 'Genome::Config::AnalysisProject::ModelBridge',
            reverse_as => 'analysis_project',
            is_many => 1,
        },
        models => {
            is => 'Genome::Model',
            via => 'model_bridges',
            to => 'model',
            is_many => 1,
        },
        config_items => {
            is => 'Genome::Config::Profile::Item',
            reverse_as => 'analysis_project',
            is_many => 1,
        },
    ],
    has_transient_optional => [
        configuration_profile => { is => 'Genome::Config::Profile' },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
};

UR::Observer->register_callback(subject_class_name => __PACKAGE__, aspect => 'is_cle', callback => \&_is_updated);
UR::Observer->register_callback(subject_class_name => __PACKAGE__, aspect => 'status', callback => \&_is_updated);

sub __display_name__ {
    my $self = shift;
    return sprintf('%s (%s)', $self->name, $self->id);
}

sub delete {
    my $self = shift;

    my $msg = 'Cannot delete analysis project %s because it has %s.';
    if ($self->model_bridge_set->count) {
        die $self->error_message($msg, $self->__display_name__, 'models');
    }
    if ($self->analysis_project_bridge_set->count) {
        die $self->error_message($msg, $self->__display_name__, 'instrument data');
    }

    my @events = Genome::Timeline::Event::AnalysisProject->get(analysis_project => $self);
    for ($self->config_items, $self->subject_mappings, @events) {
        $_->delete();
    }

    return $self->SUPER::delete();
}

sub get_configuration_profile {
    my $self = shift;

    unless ($self->configuration_profile) {
        $self->configuration_profile(
            Genome::Config::Profile->create_from_analysis_project($self)
        );
    }

    return $self->configuration_profile;
}

sub _is_updated {
    my ($self, $aspect, $old_val, $new_val) = @_;
    my $method = $aspect eq 'is_cle' ? 'cle_changed' : 'status_changed';
    Genome::Timeline::Event::AnalysisProject->$method(
        "$old_val:$new_val",
        $self,
    );
}

sub is_current {
    my $self = shift;

    my $status = $self->status;
    return if any { $_ eq $status } (qw(Completed Archived Template Deprecated));

    return 1;
}

sub environment_config_dir {
    my $self = shift;

    my $allocation = $self->disk_allocations;
    if ($allocation) {
        my $config_path = File::Spec->join($allocation->absolute_path, Genome::Config::config_subpath);
        return $allocation->absolute_path if -e $config_path;
    }

    return;
}

sub system_analysis_project {
    my $class = shift;

    my $name = Genome::Config::get('system_analysis_project_name');
    my $anp = $class->get(name => $name);

    unless ($anp) {
        $class->fatal_message(
            'No Analysis Project found for system_analysis_project_name: %s', $name
        );
    }

    return $anp;
}

1;
