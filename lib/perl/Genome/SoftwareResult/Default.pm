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
        command => {
            is => 'Command::V2',
            is_transient => 1,
            is_optional => 1,
            doc => 'the command from which this result was generated (transient)'
        },
    ],
    doc => 'saved command results'
};

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $command = $self->command;
    if ($command->can("resolve_allocation_directory")) {
        return $command->resolve_allocation_directory();
    }
    else {
        return "model_data/result-" . $self->id;
    }
}

sub resolve_allocation_disk_group_name {
    my $self = shift;
    my $command = $self->command;
    if ($command->can("resolve_allocation_disk_group_name")) {
        return $command->resolve_allocation_disk_group_name();
    }
    else {
        $ENV{GENOME_DISK_GROUP_MODELS};
    }
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

    # set the result before execute, though the result will update
    # this allows the command to link underlying results to it
    $command->result($self);

    my $saved_output_dir;
    if ($command->can('output_dir')) {
        if ($command->can('stage_output') and $command->stage_output) {
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
                my $actual_output_dir = $self->output_dir;
                #$command->debug_message("SYMLINK $actual_output_dir $saved_output_dir");
                #Genome::Sys->create_symlink($actual_output_dir, $saved_output_dir);
                $self->output_dir($saved_output_dir);
            }
            else {
                $self->_prepare_output_directory;
            }
            $self->_promote_data;
        }
        $self->_reallocate_disk_allocation;
    }
    
    $command->is_executed(1);

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
# around '_execute_body' => Genome::SoftwareResult::Default::execute_wrapper;
#

sub execute_wrapper {
    # This is a wrapper for real execute() calls.
    # All execute() methods are turned into _execute_body at class init, 
    # so this will get direct control when execute() is called. 
    my $orig = shift;
    my $command = shift;

    # handle calls as a class method
    my $was_called_as_class_method = 0;
    if (ref($command)) {
        if ($command->is_executed) {
            Carp::confess("Attempt to re-execute an already executed command.");
        }
    }
    else {
        # called as class method
        # auto-create an instance and execute it
        $command = $command->create(@_);
        return unless $command;
        $was_called_as_class_method = 1;
    }

    # handle __errors__ objects before execute
    if (my @problems = $command->__errors__) {
        for my $problem (@problems) {
            my @properties = $problem->properties;
            $command->error_message("Property " .
                                 join(',', map { "'$_'" } @properties) .
                                 ': ' . $problem->desc);
        }
        $command->delete() if $was_called_as_class_method;
        return;
    }

    my $result;
    my $meta = $command->__meta__;

    # for now shortcutting and saving a software result isn't optional.  change?
    # tentative logic for a non-shortcut scenario is here
    # TODO: create a value-based result which is not persisted but has the same API
    #  $result = $command->_execute_body(@_); # default/normal unsaved execute
    #  $command->is_executed(1);
    #  $command->result($command);
    
    $command->debug_message("Execution preceded by check for existing software result.");
    my $result_class = $command->class . '::Result';
    my %props = _copyable_properties($command, $result_class);
    delete $props{output_dir};
    unless ($result) {
        $result = $result_class->get_or_create(
            test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
            %props, 
            command => $command,
        );
    }
    if ($command->is_executed) {
        $command->debug_message("New software result saved: " . $result->__display_name__);
    }
    else {
        $command->debug_message("Existing results found: " . $result->__display_name__);
    }

    # copy properties from the result to the command outputs/changes
    # re-check the command so that output values and metrics are copied
    %props = _copyable_properties($command, $result_class);
    for my $name (keys %props) {
        my $meta = $result->__meta__->property($name);
        next if (not defined $props{$name});
        next if $meta->is_input or $meta->is_param;
        if ($meta->is_many) {
            my @old_values = sort $result->$name;
            my @new_values = sort @{$props{$name}};
            if ("@old_values" ne "@new_values") {
                Carp::confess("has-many properties which change during execute ($name) are not currently supported until the accessor is smarter! Change from [@old_values] to [@new_values]")
            }
        }
        else {
            $result->$name($props{$name});
        }
    }
    $command->result($result);
    if ($command->can('output_dir')) {
        $command->output_dir($result->output_dir);
        $command->debug_message("Output directory: " . $result->output_dir);
    }

    return $command if $was_called_as_class_method;
    return 1;
}

sub _copyable_properties {
    my ($command,$result) = @_;
    my %props;
    my $command_meta = $command->__meta__;
    my @command_properties = $command_meta->properties();
    for my $command_property (@command_properties) {
        my $name = $command_property->property_name;
        next unless ($command_property->class_name->isa('Command'));
        if ($result->can($name)) {
            if ($command_property->is_many) {
                $props{$name} = [ $command->$name ];
            }
            else {
                $props{$name} = $command->$name;
            }
        }
    }
    return %props;
}

1;

