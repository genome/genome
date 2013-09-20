package Genome::Model::Tools::Joinx::SafeVcfMerge;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;
use YAML;

class Genome::Model::Tools::Joinx::SafeVcfMerge {
    is => 'Genome::Model::Tools::Joinx::VcfMerge',
    doc => 'Just like Joinx::VcfMerge except that it can keep you from running out of file-handles',
    has_optional_input => [
        max_files_per_merge => {
            is => 'Text',
            default => 512,
            doc => 'Set this to cause no more than N vcfs to be merged during a single operation at a time',
        },
        working_directory => {
            is => 'Path',
            doc => 'Is used to store intermediate merges (if needed).',
        },
        remove_intermediate_files => {
            is => 'Boolean',
            default => 1,
            doc => 'Remove intermediate results that were generated.',
        },
    ],
};

sub execute {
    my ($self) = @_;

    if ($self->labeled_input_files) {
        die "SafeVcfMerge does not support labeled-input-files feature";
    }

    $self->status_message("Resolving command inputs");
    if($self->max_files_per_merge < 2) {
        my $max = $self->max_files_per_merge;
        $self->error_message("max_files_per_merge=$max and it doesn't make sense to have max_files_per_merge < 2.");
        Carp::croak($self->error_message());
    }

    my @inputs = $self->_resolve_inputs();
    unless(@inputs) {
        $self->error_message("No non-empty inputs found to merge");
        Carp::croak($self->error_message());
    }

    my @input_groups = $self->_resolve_input_groups(\@inputs);
    if(scalar(@input_groups) > 1) {
        $self->status_message(sprintf("Breaking your inputs into %s groups",
                scalar(@input_groups)));
        my ($working_directory) = $self->_resolve_working_directory();
        my $joinx_log_dir = $self->_resolve_joinx_log_dir($working_directory);
        my $workflow_log_dir = $self->_resolve_workflow_log_dir($working_directory);

        # generate and execute workflow.
        $self->status_message("Generating workflow to run safe-sized merges");
        my ($workflow, $workflow_inputs, $intermediate_files) =
                $self->_generate_workflow(
                        $working_directory,
                        $joinx_log_dir,
                        $workflow_log_dir,
                        \@input_groups);

        $self->_execute_workflow($workflow, $workflow_inputs, $working_directory);

        if($self->remove_intermediate_files) {
            $self->status_message("Removing intermediate files");
            $self->_remove_intermediate_files($intermediate_files);
        }
    } else {
        ## NOTE the following cannot be used because magic on execute in Command::V1 and V2
        #return $self->SUPER::execute(@_);
        $self->status_message("Running without building a workflow");
        my $sub = $self->super_can('_execute_body');
        $sub->($self, @_);
    }
    $self->status_message("Everything completed");
    return 1;
}

# Return the non-empty inputs grouped into bunches no larger than
# max_files_per_merge.
sub _resolve_input_groups {
    my ($self, $inputs) = @_;
    my @inputs = @{$inputs};

    return _chunk_array(\@inputs, $self->max_files_per_merge);
}

# Return an array of array_refs that break the input_array into chunks with
# at most num_elements_per_chunk.
#     chunk_array([1, 2, 3, 4, 5, 6], 5)
#     returns: ([1, 2, 3, 4, 5], [6])
sub _chunk_array {
    my ($input_array_ref, $num_elements_per_chunk) = @_;
    my @input_array = @{$input_array_ref};
    my $max_index = scalar(@input_array) - 1;

    my @chunks;
    my $count = 0;
    my $high = 0;
    while($high != $max_index) {
        my $low = $count * $num_elements_per_chunk;
        $high = ($count+1) * $num_elements_per_chunk - 1;
        $high = $max_index if $high > $max_index;

        my @chunk = @input_array[$low..$high];
        push(@chunks, \@chunk);
        $count += 1;
    }
    return @chunks;
}

sub _get_uuid_string {
    my $ug = Data::UUID->new();
    my $uuid = $ug->create();
    return $ug->to_string($uuid);
}

sub _resolve_working_directory {
    my ($self) = @_;

    if(defined($self->working_directory)) {
        unless(-d $self->working_directory) {
            Genome::Sys->create_directory($self->working_directory);
            unless(-d $self->working_directory) {
                $self->error_message(
                    sprintf("Working directory (%s) does not exist" .
                            " and couldn't be created.",
                        $self->working_directory));
                Carp::croak($self->error_message());
            }
        }
        return ($self->working_directory);
    } else {
        $self->error_message("You must supply a working-directory to store" .
            " intermediate results.");
        Carp::croak($self->error_message());
    }
}

sub _create_dir {
    my ($base_path, $name) = @_;

    my $new_dir = join("/", $base_path, $name);
    unless(-d $new_dir) {
        Genome::Sys->create_directory($new_dir);
    }
    return $new_dir;
}

sub _resolve_joinx_log_dir {
    my ($self, $base_path) = @_;

    if($self->error_log) {
        $self->warning_message(
            sprintf("You specified a location (%s) for joinx errors to be" .
                    " logged but since a workflow will be launched, the" .
                    " errors from each run of joinx will be placed in the" .
                    " working directory (%s)", $self->error_log, $base_path));
    }
    return _create_dir($base_path, "joinx_logs");
}

sub _resolve_workflow_log_dir {
    my ($self, $base_path) = @_;

    return _create_dir($base_path, "logs");
}

sub inherited_joinx_property_names {
    my $self = shift;
    my @vcf_merge_properties = Genome::Model::Tools::Joinx::VcfMerge->__meta__->properties;

    # this excludes things from higher up the class heirarchy (such as UR::Object)
    my @not_is_many = grep {!$_->is_many} @vcf_merge_properties;
    my @joinx_properties = grep {$_->class_name =~ m/Joinx/} @not_is_many;
    return map {$_->property_name} @joinx_properties;
}

sub _workflow_inputs {
    my $self = shift;

    my %workflow_inputs;
    for my $inherited_property_name ($self->inherited_joinx_property_names) {
        if (defined($self->$inherited_property_name) and $inherited_property_name ne 'output_file') {
            $workflow_inputs{$inherited_property_name} = $self->$inherited_property_name;
        }
    }
    return %workflow_inputs;
}

sub _generate_workflow {
    my ($self, $working_directory, $joinx_log_dir, $workflow_log_dir,
            $input_groups) = @_;
    my @input_groups = @{$input_groups};


    my %workflow_inputs = $self->_workflow_inputs;
    my %common_workflow_inputs = %workflow_inputs;
    $workflow_inputs{final_output_file} = $self->output_file;
    $workflow_inputs{final_merge_working_directory} = join("/",
            $working_directory, "sub_working_directory");
    $workflow_inputs{final_error_log} = join("/",
            $joinx_log_dir, "final_merge.err");

    my $ext = 'vcf';
    $ext .= '.gz' if $self->use_bgzip;
    my @converge_inputs;
    my @intermediate_files;
    for my $i (0..$#input_groups) {
        $workflow_inputs{"input_files_$i"} = $input_groups[$i];

        my $group_output = sprintf("%s/output_%s.%s",
                $working_directory, $i, $ext);
        push(@intermediate_files, $group_output);
        $workflow_inputs{"output_file_$i"} = $group_output;

        my $error_log = join("/", $working_directory, "joinx_logs",
                "group_$i.err");
        $workflow_inputs{"joinx_logs_$i"} = $error_log;
        push(@converge_inputs, "Group $i");
    }

    my $workflow = Workflow::Model->create(
        name => 'Vcf Merge',
        input_properties => [sort(keys %workflow_inputs)],
        output_properties => [ "merged_vcf_file"],
    );
    $workflow->log_dir($workflow_log_dir);

    my $converge_operation = $workflow->add_operation(
        name => "Converge Inputs for Final Merge",
        operation_type => Workflow::OperationType::Converge->create(
            input_properties => \@converge_inputs,
            output_properties => [qw(input_files)],
        ),
    );

    # final merge operation
    my $final_merge_operation = $workflow->add_operation(
        name => "Final Merge",
        operation_type => Workflow::OperationType::Command->get(
                $self->class),
    );

    my %properties = (
        final_output_file => "output_file",
        final_merge_working_directory => "working_directory",
        final_error_log => "error_log"
    );
    _make_links($workflow, \%properties, $final_merge_operation);
    _connect_common_inputs($workflow, \%common_workflow_inputs,
            $final_merge_operation);

    $workflow->add_link(
        left_operation => $converge_operation,
        left_property => "input_files",
        right_operation => $final_merge_operation,
        right_property => "input_files",
    );

    $workflow->add_link(
        left_operation => $final_merge_operation,
        left_property => "output_file",
        right_operation => $workflow->get_output_connector(),
        right_property => "merged_vcf_file",
    );

    # group merge operations
    my $input_connector = $workflow->get_input_connector();
    for my $i (0..$#input_groups) {
        my $group_merge_operation = $workflow->add_operation(
            name => "Merge Group $i",
            operation_type => Workflow::OperationType::Command->get(
                    $self->class),
        );

        my %properties = (
            "input_files_$i" => "input_files",
            "output_file_$i" => "output_file",
            "joinx_logs_$i" => "error_log",
        );
        _make_links($workflow, \%properties, $group_merge_operation);
        _connect_common_inputs($workflow, \%common_workflow_inputs,
                $group_merge_operation);

        $workflow->add_link(
            left_operation => $group_merge_operation,
            left_property => "output_file",
            right_operation => $converge_operation,
            right_property => $converge_inputs[$i],
        );
    }
    return $workflow, \%workflow_inputs, \@intermediate_files;
}

sub _dump_workflow {
    my ($workflow, $output_filename) = @_;

    my $xml_file = Genome::Sys->open_file_for_writing($output_filename);
    $workflow->save_to_xml(OutputFile => $xml_file);
    $xml_file->close();
}

sub _execute_workflow {
    my ($self, $workflow, $workflow_inputs, $working_directory) = @_;

    $self->status_message("Validating the vcf-merge workflow.");
    my @errors = $workflow->validate();
    if (@errors) {
        $self->error_message("Error validating vcf-merge workflow");
        $self->error_message(@errors);
        Carp::croak($self->error_message());
    }
    _dump_workflow($workflow,
            join("/", $working_directory, 'workflow.xml'));

    $self->status_message("Launching the vcf-merge workflow.");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %{$workflow_inputs});

    unless($result){
        $self->error_message("Vcf-merge workflow did not return correctly.");
        $self->error_message(join("\n",
                map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)));
        Carp::croak($self->error_message());
    }
    return $workflow_inputs->{final_output_file};
}

# connect the properties from the workflow's input_connector to 'target_operation'.
sub _make_links {
    my ($workflow, $properties, $target_operation) = @_;
    my %properties = %{$properties};

    my $input_connector = $workflow->get_input_connector();
    for my $left_property (keys %properties) {
        my $right_property = $properties{$left_property};
        $workflow->add_link(
            left_operation => $input_connector,
            left_property => $left_property,
            right_operation => $target_operation,
            right_property => $right_property,
        );
    }
}

# connect the common inputs to the target operation.
sub _connect_common_inputs {
    my ($workflow, $common_workflow_inputs, $target_operation) = @_;
    my %common_workflow_inputs = %{$common_workflow_inputs};

    my $input_connector = $workflow->get_input_connector();
    for my $key (keys %common_workflow_inputs) {
        $workflow->add_link(
            left_operation => $input_connector,
            left_property => $key,
            right_operation => $target_operation,
            right_property => $key,
        );
    }
}

sub _remove_intermediate_files {
    my ($self, $intermediate_files) = @_;
    my @intermediate_files = @{$intermediate_files};

    for my $file (@intermediate_files) {
        unlink($file);
    }
}

1;
