package Genome::Command::WithSoftwareResult;
use strict;
use warnings;


class Genome::Command::WithSoftwareResult {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        _software_result_version => {
            is_param => 1,
            is => 'Integer',
            is_abstract => 1,
            doc => 'the version of software_results generated, which may iterate as this logic iterates'
        },
    ],
    has_transient => [
        initial_allocation_size => {
            default_value => 50_000_000,
            doc => 'the initial size of the disk allocation if <stage_output> is false',
        },
        software_result => {
            is => 'Genome::SoftwareResult',
            is_optional => 1,
            doc => 'the saved software result'
        },
    ],
    is_abstract => 1,
};

# this is where your execute logic belongs now
sub _execute {
    my ($self, $output_dir) = @_;
    die "Abstract";
}


# ===== OPTIONAL STUFF TO OVERWRITE =====

# overwrite this to do something after the software result has been
# created (for instance, symlink results into place).
sub _finalize {
    my ($self, $result) = @_;
    return;
}

# return one of the keys in %SOFTWARE_RESULT_TYPES (found below)
sub software_result_type {
    return 'basic';
}

# disk allocation information
sub resolve_allocation_subdirectory {
    my ($self, $result) = @_;
    return "model_data/result-" . $result->id;
}

# disk allocation information
sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

# ==== INTERNAL USE ONLY ====
sub execute {
    my $self = shift;

    $self->status_message("Execution preempted by check for existing software result.");
    my $result_class = $self->class . '::Result';
    my %args = _params_and_inputs($self);
    $args{test_name} = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef;
    $args{command} = $self;

    $self->status_message("Searching for software result with args: " . Data::Dumper::Dumper(%args));
    my $result = $result_class->get_with_lock(%args);
    if ($result) {
        $self->status_message("Existing results found: " . $result->__display_name__);
        $self->software_result($result);
    } else {
        $self->status_message("No existing software result found, creating one.");
        $result = $result_class->create(%args);
        $self->status_message("New software result saved: " . $result->__display_name__);
    }

    $self->_finalize($result);
    return 1;
}

sub shortcut {
    my $self = shift;

    $self->status_message("Attempting to shortcut.");
    my $result_class = $self->class->result_class();
    my %args = _params_and_inputs($self);
    $args{test_name} = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef;
    $args{command} = $self;

    my $result = $result_class->get_with_lock(%args);
    if ($result) {
        $self->status_message("Existing result found: " . $result->__display_name__);
        $self->software_result($result);
        $self->_finalize($result);
        return 1;
    } else {
        $self->status_message("No existing software result found.  Shortcutting failed.");
        return;
    }
}


sub _params_and_inputs {
    my $self = shift;

    my %props;
    my $meta = $self->__meta__;
    my @properties = $meta->properties();
    for my $property (@properties) {
        my $name = $property->property_name;
        my $is_param = $property->can('is_param') && $property->is_param;
        my $is_input = $property->can('is_input') && $property->is_input;
        if ($is_param or $is_input) {
            if ($property->is_many) {
                $props{$name} = [ $self->$name ];
            }
            else {
                $props{$name} = $self->$name;
            }
        }
    }
    return %props;
}


# === all magic below this point ====
our %SOFTWARE_RESULT_TYPES = (
    staged => 'Genome::SoftwareResult::WithCommand::Staged',
    basic => 'Genome::SoftwareResult::WithCommand::Basic',
);

sub __extend_namespace__ {
    my ($self, $ext) = @_;

    my $meta = $self->SUPER::__extend_namespace__($ext);
    if ($meta) {
        return $meta;
    }

    my $new_class = $self . '::' . $ext;
    my ($command_class, $final_ext) = ($new_class =~ /^(.*)::(.*)$/);

    # ensure the command class is real
    my $cmd_meta = UR::Object::Type->get($command_class);
    return unless $cmd_meta;

    my $base_class = $SOFTWARE_RESULT_TYPES{$command_class->software_result_type};
    return unless $base_class;

    my %has = Command::V2::_wrapper_has($command_class, $base_class);

    class {$new_class} {
        is => $base_class,
        has => [ %has ],
    };
}


1;
