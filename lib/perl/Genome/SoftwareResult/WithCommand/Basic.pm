package Genome::SoftwareResult::WithCommand::Basic;
use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::WithCommand::Basic {
    is => ['Genome::SoftwareResult'],
    has => [
        command => {
            is => 'Genome::Command::WithSoftwareResult',
            is_transient => 1,
            is_abstract => 1,
            doc => 'the command from which this result was generated (transient)'
        },
    ],
    is_abstract => 1,
};

# just delegate back to the command for these
sub resolve_allocation_subdirectory {
    my $self = shift;
    return $self->command->resolve_allocation_subdirectory($self);
}


sub resolve_allocation_disk_group_name {
    my $self = shift;
    return $self->command->resolve_allocation_disk_group_name($self);
}


sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    # create disk allocation
    my $allocation = $self->_create_disk_allocation();
    my $output_dir = $allocation->absolute_path;

    my $command = $self->command;
    $command->software_result($self);
    $command->_execute($output_dir);

    $allocation->reallocate();

    return $self;
}


sub _create_disk_allocation {
    my $self = shift;

    my %allocation_params = (
        disk_group_name => $self->resolve_allocation_disk_group_name,
        allocation_path => $self->resolve_allocation_subdirectory,
        kilobytes_requested => $self->command->initial_allocation_size,
        owner_class_name => $self->class,
        owner_id => $self->id,
    );

    my $allocation = Genome::Disk::Allocation->allocate(%allocation_params);
    unless ($allocation) {
        $self->error_message("Failed to get disk allocation with params:\n". Data::Dumper::Dumper(%allocation_params));
        die($self->error_message);
    }

    my $output_dir = $allocation->absolute_path;
    $self->output_dir($output_dir);

    return $allocation;
}

1;
