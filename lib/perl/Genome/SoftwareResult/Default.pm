package Genome::SoftwareResult::Default;
use strict;
use warnings;
use Genome;

# Genome::SomeCommand::Result is auto-generated
# with this as its base class for any command 
# in this namespace.
# To get it to automatically use it set this in 
# any command class definition:
#
#   shortcut_execute => 1

class Genome::SoftwareResult::Default {
    is => 'Genome::SoftwareResult::Stageable',
    has => [
        command  => {
            is => 'Command::V2',
            is_transient => 1,
            doc => 'the command from which this result was generated (transient)'
        },
    ],
    doc => 'saved command results'
};

sub resolve_allocation_subdirectory {
    my $self = shift;
    return "model_data/result-" . $self->id;
}

sub resolve_allocation_disk_group_name {
    "info_genome_models" 
}

sub create {
    my $class = shift;

    if ($class eq __PACKAGE__ or $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    my $bx = $class->define_boolexpr(@_);
    my $command = $bx->value_for('command');
    unless ($command) {
        # not creating based on a command (backfill?)
        return $class->SUPER::create(@_);
    }

    # copy properties from the command at construction time
    # as they are essential to creating the correct lock
    my $self = $class->SUPER::create(@_);
    return unless $self;

    my $saved_output_dir;
    if ($command->can('output_dir')) {
        if ($command->__meta__->can('stage_output') and $command->__meta__->stage_output) {
            $self->_prepare_staging_directory;
            $saved_output_dir = $command->output_dir;
            $command->output_dir($self->temp_staging_directory);
        }
        elsif (not $command->output_dir) {
            $self->_prepare_output_directory;
            $command->output_dir($self->output_dir);
        }
        else {
            # output dir is as specified by the caller on the command
        }
        die "no output dir set???" unless $command->output_dir;
    }

    $command->_execute_body();

    if ($command->output_dir) {
        if ($self->temp_staging_directory) {
            if ($saved_output_dir) {
                $self->output_dir($saved_output_dir);
            }
            else {
                $self->_prepare_output_directory;
            }
            $self->_promote_data;
        }
        $self->_reallocate_disk_allocation;
    }

    return $self;
}

sub _staging_disk_usage {
    my $self = shift;

    my $tmp = $self->temp_staging_directory;
    unless ($tmp) {
        # TODO: delegate through command for a better estimate
        return 1;
    }

    my $usage;
    unless ($usage = Genome::Sys->disk_usage_for_path($self->temp_staging_directory)) {
        $self->error_message("Failed to get disk usage for staging: " . Genome::Sys->error_message);
        die $self->error_message;
    }

    return $usage;
}


1;
