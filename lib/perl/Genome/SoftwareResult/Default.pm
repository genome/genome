package Genome::SoftwareResult::Default;
use strict;
use warnings;
use Genome;

# This is the base class for auto-generated software result subclasses
# The Genome.pm __extend_namespace__ creates these with a class name
# like ${COMMANDCLASS}::Result for any Command::V2 in a consistent way.

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

#
# Add this to a Command to get automatic default software results
#
# use Moose;
# around 'execute' => Genome::SoftwareResult::Default::execute_wrapper;
#

sub execute_wrapper {
    # This is a wrapper for real execute() calls.
    # All execute() methods are turned into _execute_body at class init, 
    # so this will get direct control when execute() is called. 
    my $orig = shift;
    my $self = shift;

    # handle calls as a class method
    my $was_called_as_class_method = 0;
    if (ref($self)) {
        if ($self->is_executed) {
            Carp::confess("Attempt to re-execute an already executed command.");
        }
    }
    else {
        # called as class method
        # auto-create an instance and execute it
        $self = $self->create(@_);
        return unless $self;
        $was_called_as_class_method = 1;
    }

    # handle __errors__ objects before execute
    if (my @problems = $self->__errors__) {
        for my $problem (@problems) {
            my @properties = $problem->properties;
            $self->error_message("Property " .
                                 join(',', map { "'$_'" } @properties) .
                                 ': ' . $problem->desc);
        }
        $self->delete() if $was_called_as_class_method;
        return;
    }

    my $result;
    my $meta = $self->__meta__;
    if ($meta->shortcut_execute) {
        my $results_class = $meta->save_results_as;
        unless ($results_class) {
            die "cannot shortcut for " . $self->class . " without a setting for save_result_as!";
        }
        #unless ($results_class->can('get_or_create_for_command')) {
        #    die "class $results_class does not implement get_or_create_for_command()";
        #}
        #$result = $self->shortcut;
        unless ($result) {
            my %props = $self->_copyable_properties_for($results_class);
            $result = $results_class->get_or_create(%props, command => $self);
        }

        # copy properties from the result to the command outputs/changes
        my %props = $self->_copyable_properties_for($results_class);
        for my $name (keys %props) {
            $result->$name($props{$name});
        }
    }
    else {
        if (my $results_class = $meta->save_results_as) {
            die "save_results_as without shortcut_execute is not currently implemented!";
        }
        $result = $self->_execute_body(@_); # default/normal unsaved execute
        $self->is_executed(1);
    }
    $self->result($result);


    return $self if $was_called_as_class_method;
    return $result;
}
1;

