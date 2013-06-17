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
    my ($self, $result) = @_;
    return "info_genome_models";
}

# ==== INTERNAL USE ONLY ====
sub execute {
    my $self = shift;

    $self->status_message("Execution preempted by check for existing software result.");
    my $result_class = $self->class . '::Result';
    my %args = _params_and_inputs($self);
    $args{test_name} = $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef;
    $args{command} = $self;

    my $result = $result_class->get_with_lock(%args);
    if ($result) {
        $self->status_message("Existing results found: " . $result->__display_name__);
        $self->software_result($result);
        #_copy_outputs($result, $self);
    } else {
        $self->status_message("No existing software result found, creating one.");
        $result = $result_class->create(%args);
        #_ensure_inputs_and_params_didnt_change_during_execution(); # TODO
        $self->status_message("New software result saved: " . $result->__display_name__);
        #_copy_outputs($self, $result);
    }

    $self->_finalize($result);
    return 1;
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

# NEED TO WAIT FOR A TABLE TO BE MADE TO DO THIS STUFF PROPERLY
#sub _copy_outputs {
#    my ($source, $destination) = @_;
#
#    my %props = _copyable_properties($source, $destination);
#    for my $name (keys %props) {
#        my $meta = $source->__meta__->property($name);
#        if ($meta->can('is_output') && $meta->is_output) {
#            if ($meta->is_many) {
#                my @values = $source->$name;
#                $destination->$name(\@values);
#            } else {
#                $destination->$name($source->$name);
#            }
#        }
#    }
#}
#
#
#our %DO_NOT_COPY = (
#    id => 1,
#    _lock_name => 1,
#);
#
#sub _copyable_properties {
#    my ($source, $destination) = @_;
#    my %props;
#    my $source_meta = $source->__meta__;
#    my @source_properties = $source_meta->properties();
#    for my $source_property (@source_properties) {
#        my $name = $source_property->property_name;
#        next if exists($DO_NOT_COPY{$name});
#        if ($destination->can($name)) {
#            if ($source_property->is_many) {
#                $props{$name} = [ $source->$name ];
#            }
#            else {
#                $props{$name} = $source->$name;
#            }
#        }
#    }
#    return %props;
#}

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


sub _init_subclass {
    my $subclass_name = shift;

    my @problems;
    my $meta = $subclass_name->__meta__;

    my $result_version_meta = $meta->property("_software_result_version");
    unless ($result_version_meta) {
        push @problems, "$subclass_name should implement _software_result_version, typically with a default_value of '1'";
    }

    if (@problems) {
        for (@problems) { $subclass_name->error_message($_)  }
        die "error defining $subclass_name";
    }

    return 1;
}


1;
