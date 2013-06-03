package Genome::SoftwareResult::WithCommand::Staged;
use strict;
use warnings;
use Genome;

class Genome::SoftwareResult::WithCommand::Staged {
    is => ['Genome::SoftwareResult::Stageable'],
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

    my $command = $self->command;

    my $output_dir = $self->_prepare_staging_directory();

    $command->software_result($self);
    $command->_execute($output_dir);

    $self->_prepare_output_directory();
    $self->_promote_data();
    $self->_reallocate_disk_allocation();

    return $self;
}


1;
