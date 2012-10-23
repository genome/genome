package Genome::DataSource::Main;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::DataSource::Main {
    is => 'UR::DataSource::Pg',
    has_constant => [
        server => { default_value => 'dbname=genome' },
        login => { default_value => 'genome' },
        auth => { default_value => undef },
        owner => { default_value => 'public' },
    ],
};

sub _ds_tag {
    'Genome::DataSource::Main';
}

sub _sync_database {
    my $self = shift;
    my %params = @_;

    # Need to remove all commit observers, they will be fired during commit and that's no good!
    my @observers = UR::Observer->get();
    for my $observer (@observers) {
        $observer->delete;
    }

    # Need to update all classes that have been changed (and all of their parent
    # classes) to use the postgres datasource instead of Oracle.
    my %classes = map { $_->class => 1 } @{$params{changed_objects}};
    for my $class (sort keys %classes) {
        my $meta = UR::Object::Type->get($class);
        next unless $meta;
        my @metas = ($meta, $meta->ancestry_class_metas);
        for my $meta (@metas) {

            # If the object has been deleted and we're dealing with the Ghost, need to find the non-Ghost class
            # and add that to the list of metas we need to update. Otherwise, that class won't get updated to use
            # postgres and an error will occur.
            if ($meta->class_name =~ /::Ghost$/) {
                my $non_ghost_class = $meta->class_name;
                $non_ghost_class =~ s/::Ghost$//;
                my $non_ghost_meta = UR::Object::Type->get($non_ghost_class);
                if ($non_ghost_meta) {
                    push @metas, $non_ghost_meta;
                }
                else {
                    Carp::confess "Could not find meta object for non-ghost class $non_ghost_class!";
                }
            }

            if (defined $meta->data_source) {
                $self->rewrite_classdef_to_use_postgres($meta);     
            }
        }
    }

    # Update meta datasource to point to an empty file so we don't get failures due to not being
    # able to find column/table meta data.
    my $temp_file_fh = File::Temp->new;
    my $temp_file = $temp_file_fh->filename;
    $temp_file_fh->close;

    my $meta_ds = Genome::DataSource::Meta->_singleton_object;
    $meta_ds->get_default_handle->disconnect;
    $meta_ds->server($temp_file);
    
    return $self->SUPER::_sync_database(@_);
}

# called whenever we generate an ID
sub autogenerate_new_object_id_for_class_name_and_rule {
    my $self = shift;
    UR::Object::Type->autogenerate_new_object_id_uuid; 
}

# called before attempting to write an SQL query
our %rewritten;
sub _generate_class_data_for_loading {
    my ($self, $meta) = @_;
    my @ancestor_metas = $meta->ancestry_class_metas;
    for my $m ($meta, @ancestor_metas) {
        $self->rewrite_classdef_to_use_postgres($m);
    }
    $self->SUPER::_generate_class_data_for_loading($meta);
}

# used before query _and_ from _sync_database
sub rewrite_classdef_to_use_postgres {
    my $self = shift;
    my $meta = shift;

    if ($rewritten{$meta->id}) {
        return;
    }
    $rewritten{$meta->id} = 1;
    
    my $class = $meta->class_name;

    $meta->data_source($self->_ds_tag);

    # Columns are stored directly on the meta object as an optimization, need to be updated
    # in addition to table/column objects.
    my (undef, $cols) = @{$meta->{'_all_properties_columns'}};
    $_ = lc $_ foreach (@$cols);

    if (defined $meta->table_name) {
        my $oracle_table = $meta->table_name;
        my $postgres_table = $self->postgres_table_name_for_oracle_table($oracle_table);
        unless ($postgres_table) {
            Carp::confess "Could not find postgres equivalent for oracle table $oracle_table while working on class $class";
        }
        $meta->table_name($postgres_table);

        my @properties = $meta->all_property_metas;
        for my $property (@properties) {
            next unless $property->column_name;
            $property->column_name(lc $property->column_name);
        }
    }
}

sub postgres_table_name_for_oracle_table {
    my $self = shift;
    my $oracle_table = shift;
    return unless $oracle_table;
    my %mapping = $self->oracle_to_postgres_table_mapping;
    return $mapping{lc $oracle_table};
}

sub oracle_to_postgres_table_mapping {
    return (
        'search_index_queue' => 'web.search_index_queue',
        'feature_list' => 'model.feature_list',
        'fragment_library' => 'instrument.fragment_library',
        'genome_disk_allocation' => 'disk.allocation',
        'disk_volume' => 'disk.volume',
        'disk_group' => 'disk.group',
        'disk_volume_group' => 'disk.volume_group_bridge',
        'genome_model' => 'model.model',
        'genome_model_build' => 'model.build',
        'genome_model_build_input' => 'model.build_input',
        'genome_model_build_link' => 'model.build_link',
        'genome_model_metric' => 'model.build_metric',
        'genome_model_event' => 'model.event',
        'genome_model_event_input' => 'model.event_input',
        'genome_model_event_metric' => 'model.event_metric',
        'genome_model_event_output' => 'model.event_output',
        'genome_model_group' => 'model.model_group_bridge',
        'genome_model_input' => 'model.model_input',
        'genome_model_link' => 'model.model_link',
        'genome_nomenclature' => 'web.nomenclature',
        'genome_nomenclature_enum_value' => 'web.nomenclature_enum_value',
        'genome_nomenclature_field' => 'web.nomenclature_field',
        'genome_project' => 'subject.project',
        'genome_project_part' => 'subject.project_part',
        'genome_subject' => 'subject.subject',
        'genome_subject_attribute' => 'subject.subject_attribute',
        'genome_sys_user' => 'subject.user',
        'genome_sys_user_role' => 'subject.user_role',
        'genome_sys_user_role_member' => 'subject.user_role_member',
        'genome_task' => 'web.task',
        'genome_task_params' => 'web.task_params',
        'instrument_data' => 'instrument.data',
        'instrument_data_attribute' => 'instrument.data_attribute',
        'misc_attribute' => 'subject.misc_attribute',
        'misc_note' => 'subject.misc_note',
        'model_group' => 'model.model_group',
        'processing_profile' => 'model.processing_profile',
        'processing_profile_param' => 'model.processing_profile_param',
        'software_result' => 'result.software_result',
        'software_result_input' => 'result.input',
        'software_result_metric' => 'result.metric',
        'software_result_param' => 'result.param',
        'software_result_user' => 'result.user',
        'genome_model_build_variant' => 'model.build_variant',
        'genome_model_variant' => 'model.variant',
    );
}

sub commit {
    my $self = shift;
    return $self->SUPER::commit(@_);
}
1;

