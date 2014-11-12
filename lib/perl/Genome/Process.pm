package Genome::Process;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;
use Params::Validate qw(validate_pos :types);
use Data::Dump qw(pp);

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
    ],
    doc => 'A base class to manage meta-data related to running a process (workflow)',
};


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    $self->update_status('New');
    $self->software_revision(Genome::Sys->snapshot_revision);

    return $self;
}

my $SET_TIMESTAMP_ON_STATUS = {
    New => 'created_at',
    Running => 'started_at',
    Crashed => 'ended_at',
    Succeeded => 'ended_at',
};

sub update_status {
    my ($self, $new_status) = validate_pos(@_, OBJECT, SCALAR);

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
        $self->status($new_status);
    } else {
        die sprintf("Cannot transition Process (%s) from (%s) to (%s)",
            $self->id, $old_status, $new_status);
    }
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


1;
