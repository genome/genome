package Genome::Process;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;
use Params::Validate qw(validate_pos :types);
use Data::Dump qw(pp);
use Scalar::Util qw();

class Genome::Process {
    is => [
        "Genome::Utility::ObjectWithCreatedBy",
    ],
    is_abstract => 1,
    table_name => 'process.process',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',

    subclassify_by => 'subclass_name',
    id_generator => '-uuid',
    id_by => [
        id => {is => 'Text'},
    ],

    subclass_description_preprocessor => 'Genome::Process::_preprocess_for_inputs',
    attributes_have => [
        is_input => {
            is => 'Boolean',
            default => 0,
        }
    ],

    has => [
        disk_allocation => {
            is => 'Genome::Disk::Allocation',
            id_by => 'disk_allocation_id',
            doc => 'A disk allocation to store logs and other process meta-data',
        },
        status => {
            is => 'Text',
            doc => 'The current status of the process',
        },
        status_events => {
            is => 'Genome::Process::StatusEvent',
            where => [ -order_by => [ "timestamp" ] ],
            is_many => 1,
            reverse_as => 'process',
        },
        created_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'New'",
        },
        started_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'Running'",
        },
        ended_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'Crashed' or 'Succeeded'",
        },
        software_revision => {
            is => 'Text',
            is_mutable => 0,
            doc => 'The version of the Genome code the process was created with',
        },
        results => {
            is => 'Genome::SoftwareResult',
            via => 'result_users',
            to => 'software_result',
        },
        result_users => {
            is => 'Genome::SoftwareResult::User',
            is_many => 1,
            reverse_as => 'user',
        },
        subclass_name => {
            is => 'Text',
            is_mutable => 0,
            doc => 'Used by UR to instantiate the correct sub-class',
        },
        inputs => {
            is => 'Genome::Process::Input',
            is_many => 1,
            is_optional => 1,
            reverse_as => 'process',
        }
    ],
    doc => 'A base class to manage meta-data related to running a process (workflow)',
};


sub create {
    my $class = shift;
    # we do not support creating from a boolexpr because the boolexpr logic
    # will silently filter out some invalid param-values.
    my %params = @_;

    unless (exists $params{software_revision}) {
        $params{software_revision} = Genome::Sys->snapshot_revision();
    }

    my $self = $class->SUPER::create(
        $class->_preprocess_params_for_create(%params));
    return unless $self;

    $self->status('New');

    return $self;
}

my $SET_TIMESTAMP_ON_STATUS = {
    New => 'created_at',
    Running => 'started_at',
    Crashed => 'ended_at',
    Succeeded => 'ended_at',
};

sub status {
    my ($self, $new_status) = validate_pos(@_, OBJECT,
        {type => SCALAR, optional=> 1});
    return $self->__status unless $new_status;

    my $old_status = $self->status;

    my $now = UR::Context->current->now;
    my $event = Genome::Process::StatusEvent->create(
        process => $self,
        old_status => $old_status,
        new_status => $new_status,
        timestamp => $now,
    );

    if ($event) {
        my $timestamp_accessor = $SET_TIMESTAMP_ON_STATUS->{$new_status};
        if ($timestamp_accessor) {
            $self->$timestamp_accessor($now);
        }
        $self->__status($new_status);
    } else {
        die sprintf("Cannot transition Process (%s) from (%s) to (%s)",
            $self->id, $old_status, $new_status);
    }
    return $self->__status;
}


sub workflow_name {
    my $self = shift;
    return sprintf('Genome::Process(%s)', $self->id);
}

sub newest_workflow_instance {
    my $self = shift;
    my @sorted = sort {$b->id <=> $a->id} $self->_workflow_instances;
    if (@sorted) {
        return $sorted[0];
    } else {
        return;
    }
}

sub _workflow_instances {
    my $self = shift;
    my @instances = Workflow::Operation::Instance->get(
        name => $self->workflow_name,
    );
    return @instances;
}

sub log_directory {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'logs');
}

sub metadata_directory {
    my $self = shift;
    return $self->disk_allocation->absolute_path;
}

sub create_disk_allocation {
    my ($self, $kb_requested) = validate_pos(@_, OBJECT,
        {default => 50 * 1024},
    );

    if ($self->disk_allocation) {
        $self->status_message("Process (%s) already has a disk_allocation (%s)",
            $self->id, $self->disk_allocation->id);
        return $self->disk_allocation;
    }

    my $allocation_path = File::Spec->join('model_data', 'process', $self->id);

    my %args = (
        disk_group_name => $ENV{GENOME_DISK_GROUP_MODELS},
        allocation_path => $allocation_path,
        kilobytes_requested => $kb_requested,
        owner_class_name => $self->class,
        owner_id => $self->id,
    );
    $self->debug_message("Attempting to create a disk_allocation for process " .
        "(%s) with arguments: %s", $self->id, pp(\%args));
    my $disk_allocation = Genome::Disk::Allocation->create(%args);

    if ($disk_allocation) {
        $self->disk_allocation($disk_allocation);
        $self->status_message("Process (%s) now has disk_allocation (%s)",
            $self->id, $self->disk_allocation->id);
        return $self->disk_allocation;
    } else {
        die sprintf("Failed to create disk allocation with " .
            "arguments: %s", pp(\%args));
    }
}


sub _preprocess_for_inputs {
    my ($class, $class_description) = @_;

    while (my ($property_name, $property_description) =
            each %{$class_description->{has}}) {
        if ($property_description->{is_input}) {
            _expand_is_input($class, $property_description, $property_name);
        }
    }
    return $class_description;
}

sub _expand_is_input {
    my $class = shift;
    my $property_description = shift;
    my $property_name = shift;

    $property_description->{'is_mutable'} = 0;
    $property_description->{'via'} = 'inputs';
    $property_description->{'is_delegated'} = 1;
    $property_description->{'to'} = 'value_id';
    $property_description->{'where'} = [
        name => $property_name,
        '-order_by' => 'array_index',
    ];
    if (my $data_type = $property_description->{data_type}) {
        my $value_class_name = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
            $data_type, $class);
        unless ($value_class_name->isa("UR::Value")) {
            $property_description->{'to'} = 'value';
        }
    }
}

sub _preprocess_params_for_create {
    my $self = shift;
    my %params = @_;

    my $input_params_list = $params{inputs} || [];
    for my $property ($self->__meta__->properties) {
        my $property_name = $property->property_name;

        if ($property->{is_input}) {
            my $value = delete $params{$property_name};
            next unless defined($value);
            if ($property->{is_many}) {
                push @$input_params_list,
                    @{_get_input_params_list($property, $value)};
            } else {
                push @$input_params_list, _get_input_params($property, $value);
            }
        }
    }
    $params{inputs} = $input_params_list;
    return %params;
}

sub _get_input_params_list {
    my ($property, $value_list) = @_;

    if (!ref $value_list || (ref $value_list && ref $value_list ne 'ARRAY')) {
        die sprintf(
            "Cannot set 'is_many' input named (%s) with non-arrayref (%s)",
            $property->property_name, $value_list);
    }

    my $input_params_list = [];
    if (defined($value_list) && scalar(@$value_list)) {
        my $index = 0;
        for my $value (@$value_list) {
            push @$input_params_list,
                _get_input_params($property, $value, $index);
            $index++;
        }
    }
    return $input_params_list;
}

sub _get_input_params {
    my ($property, $value, $index) = validate_pos(@_, 1, 1, 0);

    my $input_params = {name => $property->property_name};
    $input_params->{array_index} = $index if defined($index);

    if ($property->to eq 'value') {
        if (Scalar::Util::blessed($value) &&
                $property->is_valid_storage_for_value($value)) {
            $input_params->{value} = $value;
        } else {
            die sprintf("Cannot set input named (%s) to value (%s), " .
                "because it is not a (%s)",
                $property->property_name, $value, $property->data_type);
        }
    } else {
        if (ref($value)) {
            die sprintf("Cannot set input named (%s) to value (%s), " .
                "because it is a reference (%s)",
                $property->property_name, $value, ref($value));
        } else {
            $input_params->{value_class_name} = 'UR::Value',
            $input_params->{value_id} = $value,
        }
    }
    return $input_params;
}


1;
