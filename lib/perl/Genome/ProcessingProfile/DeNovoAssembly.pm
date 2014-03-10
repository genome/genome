package Genome::ProcessingProfile::DeNovoAssembly;

use strict;
use warnings;

use Genome;

use Genome::Model::DeNovoAssembly::SxReadProcessor;
use Regexp::Common;
use YAML;

class Genome::ProcessingProfile::DeNovoAssembly {
    is => 'Genome::ProcessingProfile',
    has_param => [
        coverage => {
            is => 'Number',
            is_optional => 1,
            doc => 'Use genome size to limit the number of reads used in the assembly to obtain this coverage.',
        },
        # Assembler
        assembler_name => {
            doc => 'Name of the assembler.',
            valid_values => [
                # Normal assemblies
                'abyss parallel',
                'allpaths de-novo-assemble',
                'newbler de-novo-assemble',
                'soap de-novo-assemble',
                'velvet one-button',

                # Imports
                'allpaths import',
                'newbler import',
                'soap dacc-download',
                'soap import',
                'velvet import'],
        },
        assembler_version => {
            doc => 'Version of assembler.',
            #dacc for soap import
        },
        assembler_params => {
            is_optional => 1,
            doc => 'A string of parameters to pass to the assembler.',
        },
        # Read Coverage, Trim and Filter
        read_processor => {
            is_optional => 1,
            doc => "String of read trimmers, filters and sorters to use. Find processors in 'gmt sx.' List each porocessor in order of execution as they would be run on the command line. Do not include 'gmt sx', as this is assumed. List params starting w/ a dash (-), followed by the value. Separate processors by a pipe w/ a space on each side ( | ). The read processors will be validated. Ex:\n\ttrim bwa-style --trim-qual-length | filter by-length filter-length 70",
        },
        #post assemble tools to run
        post_assemble => {
            is_optional => 1,
            doc => 'String of things to run in post assembly stage .. by default already run WU post assemble process .. more later',
        },
    ],
};

sub is_imported {
    my $self = shift;
    if ( $self->assembler_name =~ /download/ or $self->assembler_name =~ /import/ ) {
        return 1;
    }
    return;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless $self;

    # Read coverage
    if ( defined $self->coverage ) {
        # Gotta be an int, gt 0 and even
        unless ( $self->coverage =~ /^$RE{num}{real}$/
                and $self->coverage > 0
                and $self->coverage <= 500
        ) {
            # TODO pick a better number??
            $self->error_message(
                "Invalid coverage (".$self->coverage."). Coverage must be an integer greater than 0 and less than  501."
            );
            $self->delete;
            return;
        }
    }

    # Validate assembler & params
    unless ( $self->_validate_assembler_and_params ) {
        $self->status_message("Create failed - could not validate assembler and params");
        $self->delete;
        return;
    }

    # Validate read processor
    unless ( $self->_validate_read_processor ) {
        $self->status_message("Create failed - could not validate read processor");
        $self->delete;
        return;
    }

    #validate post assemble steps
    if ( $self->post_assemble ) {
        unless ( $self->_validate_post_assemble_steps ) {
            $self->status_message("Failed to validate post assemble steps");
            $self->delete;
            return;
        }
    }

    return $self;
}

sub assembler_accessor_name { #returns soap_de_novo_assemble from 'soap de-novo-assemble'
    my $self = shift;

    my $name = $self->assembler_name;
    $name =~ s/ |-/_/g;

    return $name;
}

sub assembler_class {
    my $self = shift;

    my $assembler_name = $self->assembler_name;
    my ($base, $subclass) = split (/\s+/, $assembler_name);

    $subclass =~ s/-/ /g;
    $subclass = Genome::Utility::Text::string_to_camel_case( $subclass );

    my $name = 'Genome::Model::Tools::'. ucfirst $base .'::'. $subclass;

    #TODO check makes sure it exists here??

    return $name;
}

sub tools_base_class {
    my $self = shift;
    my $base_name = $self->assembler_base_name;
    return 'Genome::Model::Tools::' . ucfirst $base_name;
}

sub assembler_base_name {
    my $self = shift;
    my @tmp = split(' ', $self->assembler_name);
    return $tmp[0];
}

sub assembler_params_as_hash {
    my $self = shift;

    #assembler params specified in pp
    my $params_string = $self->assembler_params;
    return unless $params_string; # ok

    my %params = Genome::Utility::Text::param_string_to_hash($params_string);
    unless ( %params ) { # not
        Carp::confess(
            $self->error_message("Malformed assembler params: $params_string")
        );
    }

    return %params;
}

sub _validate_assembler_and_params {
    my $self = shift;

    $self->status_message("Validating assembler and params...");

    my $assembler_accessor_name = $self->assembler_accessor_name;

    my $assembler_class = $self->assembler_class;

    my %assembler_params;
    $assembler_params{version} = $self->assembler_version;

    #below params are needed for assembly but must be derived/calculated from instrument data
    #at the time of build so fake values are plugged in here to get eval to work

    my $add_param_method = $assembler_accessor_name.'_fake_params_for_eval';

    if ( $self->can( $add_param_method ) ) {
        my %fake_addl_params = $self->$add_param_method;
        #adds ins_length to 'velvet one-button' params
        %assembler_params = ( %assembler_params, %fake_addl_params );
    }

    my $clean_up_param_method = $assembler_accessor_name.'_clean_up_params_for_eval';
    if ( $self->can( $clean_up_param_method ) ) {
        #removes insert_size params from 'soap de-novo-assemble' params
        %assembler_params = $self->$clean_up_param_method( %assembler_params );
    }

    my $assembler = eval{ $assembler_class->create(%assembler_params); };
    unless ( $assembler ) {
        if ($self->assembler_params) {
            $self->error_message("Failed to construct assembler for params: "
                . $self->assembler_params);
        }
        $self->error_message($@);
        return;
    }

    $assembler->delete;

    $self->status_message("Assembler and params OK");

    return 1;
}

#< Read Processor >#
sub _validate_read_processor {
    my $self = shift;
    $self->status_message("Validate read processor...");

    my $read_processor = $self->read_processor;
    unless ( defined $read_processor ) { # ok
        $self->status_message("No read processor to validate, skipping...");
        return 1;
    }

    $self->status_message('Read processor: '.$read_processor);
    my $sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(processor => $read_processor);
    if ( not $sx_processor ) {
        $self->error_message('Failed to validate read processor!');
        return;
    }

    $self->status_message("Read processor...OK");
    return 1;
}

sub process_instrument_data_can_parallelize {
    my $self = shift;

    my $assembler_name = $self->assembler_name;
    for my $assember_can_parallelize ( 'allpaths de-novo-assemble', 'soap de-novo-assemble' ) {
        return 1 if $self->assembler_name eq $assember_can_parallelize;
    }

    return ;
}

sub assemble_objects {
    return 1;
}

#< post assemble steps >#
sub _validate_post_assemble_steps {
    my $self = shift;

    $self->status_message("Validating post assemble steps");

    foreach my $post_assemble_part ( $self->post_assemble_parts ) {
    my ($tool_name) = $post_assemble_part =~ /^(\S+)/;
    my ($param_string) = $post_assemble_part =~ /\S+\s+(.*)/;

    $tool_name =~ s/-/ /g;

    my $class_name = Genome::Utility::Text::string_to_camel_case( $tool_name );
    my $base_class = $self->tools_base_class; #return G:M:T:Velvet, Soap, etc

    my $class = $base_class . '::' . $class_name;

    my $class_meta;
    eval { $class_meta = $class->__meta__; };
    unless ( $class_meta ) {
        $self->error_message("Can't validate tool: $class_name, this tool does not exist: $class");
        return;
    }

    my %params;
    if ( $param_string ) {
        %params = Genome::Utility::Text::param_string_to_hash( $param_string );
    }

    unless ( $self->validate_post_assemble_class_params( $class_meta, %params ) ) {
        $self->error_message("Failed to validate params for class: $class");
        return;
    }
    }

    $self->status_message("Validated post assemble steps");

    return 1;
}

sub validate_post_assemble_class_params {
    my ($self, $class_meta, %params) = @_;

    foreach my $key ( keys %params ) {
        my $property = $class_meta->property_meta_for_name( $key );
        unless ( $property ) {
            $self->error_message("Failed to validate param, $key in class, " . $class_meta->class_name);
            return;
        }

        my $value = $params{$key};
        #check value against list of valid values
        if ( $property->valid_values ) {
            unless ( grep (/^$value$/, @{$property->valid_values}) ) {
        my $valid_values;
        for my $valid_value ( @{$property->valid_values} ) {
            $valid_values .= "$valid_value ";
        }
                $self->error_message("Failed to find param $key value $value in the list of valid values: $valid_values");
                return;
            }
        }
        if ( $property->data_type eq 'Integer' ) {
            unless ( $value =~ /^$RE{num}{int}$/ ) {
                $self->error_message("Expected property data type of Integer for param, $key but got $value");
                return;
            }
        }
        elsif ( $property->data_type eq 'Boolean' ) {
            unless ( $value == 1 or $value == 0 ) { #not sure if this is the best way to check
                $self->error_message("Expected property data type of Boolean for param, $key, but got $value");
                return;
            }
        }
        elsif ( $property->data_type eq 'Number' ) {
            unless ( $value =~ /^$RE{num}{real}$/ ) {
                $self->error_message("Expected property data type of Number of param, $key, but got $value");
                return;
            }
        }
        #else is text or string.. need to check
    }

    #check for missing required param
    for my $property ( $class_meta->_legacy_properties ) {
        if ( not $property->is_optional ) {

            my $property_name = $property->property_name;

        #exceptions .. since we don't know at this point where the assembly will end up
        if ( $property_name eq 'assembly_directory' ) {
        $self->status_message("assembly_directory is a required param for tool,". $class_meta->class_name .", it will be assigned build data_directory");
        next;
        }
        #if other required param is missing, quit
            if ( not exists $params{$property_name} ) {
                $self->error_message("Failed to get required param: $property_name for class, ".$class_meta->class_name);
                return;
            }
        }
    }

    return 1
}

sub post_assemble_parts {
    my $self = shift;

    my @post_assemble_parts = split (/\;\s+|\;/, $self->post_assemble);

    unless ( @post_assemble_parts ) {
    $self->error_message("Could not find any parts to run in string: ".$self->post_assemble);
    $self->delete;
    }

    return @post_assemble_parts;
}

# Number of cpus we are allowed to use
sub get_number_of_cpus {
    my $self = shift;

    return 1 if not defined $ENV{LSB_MCPU_HOSTS};

    my @tokens = split(/\s/, $ENV{LSB_MCPU_HOSTS});
    my $cpus = 0;
    if ( not @tokens ) {
        $self->error_message('Could not split LSB_MCPU_HOSTS: '.$ENV{LSB_MCPU_HOSTS});
        return;
    }

    for ( my $i = 1; $i <= @tokens; $i += 2 ) {
        if ( $tokens[$i] !~ /^$RE{num}{int}$/ ) {
            $self->error_message('Error parsing LSB_MCPU_HOSTS ('.$ENV{LSB_MCPU_HOSTS}.'), number of cpus is not an int: '.$tokens[$i]);
            return;
        }
        $cpus += $tokens[$i];
    }

    if ( $cpus == 0 ) {
        $self->error_message('Could not get the number of cpus from LSB_MCPU_HOSTS: '.$ENV{LSB_MCPU_HOSTS});
        return;
    }

    return $cpus;
}


sub map_workflow_inputs {
    my $self = shift;

    if ($self->is_imported()) {
        return $self->_map_workflow_inputs_for_normal_import(@_);
    } else {
        return $self->_map_workflow_inputs_for_normal_assembly(@_);
    }
}

sub _map_workflow_inputs_for_normal_import {
    my $self = shift;
    my $build = shift;

    my %inputs;

    $inputs{'build'} = $build;

    return %inputs;
}

sub _map_workflow_inputs_for_normal_assembly {
    my $self = shift;
    my $build = shift;

    my %inputs;

    $inputs{'build'} = $build;
    my @instrument_data = $build->instrument_data;
    $inputs{'instrument_data'} = \@instrument_data;

    return %inputs;
}


sub _resolve_workflow_for_build {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;

    $lsf_queue ||= $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};
    $lsf_project ||= 'build' . $build->id;

    if ($self->is_imported())  {
        return $self->_resolve_workflow_for_import(
            $build, $lsf_queue, $lsf_project);
    } else {
        return $self->_resolve_workflow_for_normal_assembly(
            $build, $lsf_queue, $lsf_project);
    }
}

sub _resolve_workflow_for_import {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => ['build'],
        output_properties => ['build'],
        log_dir => $build->log_directory);

    my $input_connector = $workflow->get_input_connector();
    my $output_connector = $workflow->get_output_connector();

    my $import_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Command::Import', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});

    $workflow->add_link(
        left_operation => $input_connector, left_property => 'build',
        right_operation => $import_op, right_property => 'build');

    if ($self->post_assemble) {
        my $post_assemble_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Command::PostAssemble', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});

        $workflow->add_link(
            left_operation => $import_op, left_property => 'build',
            right_operation => $post_assemble_op, right_property => 'build');
        $workflow->add_link(
            left_operation => $post_assemble_op, left_property => 'build',
            right_operation => $output_connector, right_property => 'build');
    } else {
        $workflow->add_link(
            left_operation => $import_op, left_property => 'build',
            right_operation => $output_connector, right_property => 'build');
    }

    return $workflow;
}


sub _resolve_workflow_for_normal_assembly {
    my ($self, $build, $lsf_queue, $lsf_project) = @_;

    # stages:
        # process instrument data
            # if parallelize: use ProcessID
            # else: use PrepareID
        # assemble
            # if parallelize: MergeAndLinkSxResults
            # Assemble
            # if post_assemble: PostAssemble
            # Report

    my $workflow = Workflow::Model->create(
        name => $build->workflow_name,
        input_properties => ['build', 'instrument_data'],
        output_properties => ['report_directory'],
        log_dir => $build->log_directory);

    my $input_connector = $workflow->get_input_connector();
    my $output_connector = $workflow->get_output_connector();

    my $assemble_op = $self->_add_assembler($workflow, $build, $lsf_queue,
        $lsf_project);

    my $id_op;
    if ($self->process_instrument_data_can_parallelize) {
        $id_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Build::ProcessInstrumentData', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});
        $workflow->add_link(
            left_operation => $input_connector, left_property => 'instrument_data',
            right_operation => $id_op, right_property => 'instrument_data');

        $id_op->parallel_by('instrument_data');

        my $merge_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Build::MergeAndLinkSxResults', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});
        $workflow->add_link(
            left_operation => $id_op, left_property => 'build',
            right_operation => $merge_op, right_property => 'build');

        $workflow->add_link(
            left_operation => $merge_op, left_property => 'output_build',
            right_operation => $assemble_op, right_property => 'build');
        $workflow->add_link(
            left_operation => $merge_op, left_property => 'sx_results',
            right_operation => $assemble_op, right_property => 'sx_results');
    } else {
        $id_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Command::PrepareInstrumentData', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});

        $workflow->add_link(
            left_operation => $id_op, left_property => 'build',
            right_operation => $assemble_op, right_property => 'build');
    }

    $workflow->add_link(
        left_operation => $input_connector, left_property => 'build',
        right_operation => $id_op, right_property => 'build');


    my $report_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Command::Report', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});
    if ($self->post_assemble) {
        my $post_assemble_op = _add_operation($workflow, 'Genome::Model::DeNovoAssembly::Command::PostAssemble', {
            lsf_queue => $lsf_queue, lsf_project => $lsf_project});

        $workflow->add_link(
            left_operation => $assemble_op, left_property => 'build',
            right_operation => $post_assemble_op, right_property => 'build');
        $workflow->add_link(
            left_operation => $post_assemble_op, left_property => 'build',
            right_operation => $report_op, right_property => 'build');
    } else {
        $workflow->add_link(
            left_operation => $assemble_op, left_property => 'build',
            right_operation => $report_op, right_property => 'build');
    }

    $workflow->add_link(
        left_operation => $report_op, left_property => 'report_directory',
        right_operation => $output_connector, right_property => 'report_directory');

    return $workflow;
}

sub _add_assembler {
    my ($self, $workflow, $build, $default_lsf_queue, $lsf_project) = @_;

    my $lsf_resource = $build->resolve_assemble_lsf_resource();
    my $assemble_lsf_queue = $build->resolve_assemble_lsf_queue || $default_lsf_queue;

    my $assembler_class = 'Genome::Model::DeNovoAssembly::';
    $assembler_class .= ( $self->process_instrument_data_can_parallelize ? 'Build' : 'Command' );
    $assembler_class .= '::Assemble';

    return _add_operation($workflow, $assembler_class, {
            lsf_queue => $assemble_lsf_queue,
            lsf_project => $lsf_project,
            lsf_resource => $lsf_resource});
}

sub _add_operation {
    my ($workflow, $command_class_name, $options) = @_;

    my $name = $command_class_name;
    $name =~ s/Genome::Model::DeNovoAssembly::(Build|Command):://;

    my $operation_type = Workflow::OperationType::Command->create(
        command_class_name => $command_class_name);

    for my $key (keys %{$options}) {
        unless ($operation_type->can($key)) {
            die "Illegal parameter specified for Workflow::OperationType::Command";
        }
        my $value = $options->{$key};
        if ($value) {
            $operation_type->$key($value);
        }
    }

    unless (defined $operation_type) {
        die "Workflow::OperationType undefined for params: "
            . "command_class_name = '$command_class_name', "
            . "options -> " . YAML::Dump($options);
    }

    return $workflow->add_operation(
        name => $name, operation_type => $operation_type);
}

1;
