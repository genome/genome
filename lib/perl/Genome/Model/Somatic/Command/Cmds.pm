package Genome::Model::Somatic::Command::Cmds;

use strict;
use warnings;

use Workflow;
use Genome;

class Genome::Model::Somatic::Command::Cmds {
    is => ['Genome::Command::Base', 'Workflow::Operation::Command'],
    has_many => [
        models => { 
            shell_args_position => 1,
            is => 'Genome::Model::Somatic',
        }, 
    ],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); }
};

sub help_brief {
    "Runs the cmds pipeline for sequencing data."
}

sub help_synopsis{
    my $self = shift;
    return <<"EOS"
genome model somatic cmds --models "123 456 789" --data-directory /someplace/for/output
EOS
}

sub help_detail {
    my $self = shift;
    return <<"EOS"
This tool runs the cmds pipeline for sequencing data 
The only parameters that should be provided are --data-directory and either --model-group-id or --model-ids.
EOS
}

# We must explicitly delegate to the Workflow execute, because of multiple inheritance
sub execute {
    my $self = shift;

    return $self->Workflow::Operation::Command::_execute_body(@_);
}

sub pre_execute {
    my $self = shift;

    # this is necessary to pass into the workflow
    my @models = $self->models;
    $self->input_models(\@models);

    my %default_filenames = $self->default_filenames;
    for my $param (keys %default_filenames) {
        # set a default param if one has not been specified
        my $default_filename = $self->data_directory . "/" . $default_filenames{$param};
        unless ($self->$param) {
            $self->status_message("Param $param was not provided... generated $default_filename as a default");
            $self->$param($default_filename);
        }
    }

    # Create directories that do not already exist
    for my $dir ($self->data_directory, $self->compile_cna_output_dir, $self->merge_output_dir, $self->region_output_dir) {
        unless (-d $dir) {
            $self->status_message("$dir does not exist... creating it.");
            Genome::Sys->create_directory($dir);
        }
    }


    return 1;
}

sub default_filenames{
    my $self = shift;
   
    my %default_filenames = (
        compile_cna_output_dir => "/compiled_cna_output/",
        merge_output_dir => "/merged_output/",
        region_output_dir => "/individual_region_output/",
        table_output => "table.out",
    );

    return %default_filenames;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="CMDS sequencing data pipeline" logDir="/gsc/var/log/genome/cmds_sequencing_data_pipeline">

  <link fromOperation="input connector" fromProperty="input_models" toOperation="Compile Cna Output" toProperty="models" />
  <link fromOperation="input connector" fromProperty="compile_cna_output_dir" toOperation="Compile Cna Output" toProperty="output_dir" />
  
  <link fromOperation="input connector" fromProperty="merge_output_dir" toOperation="Merge Cna Output By Chrom" toProperty="output_dir" />
  <link fromOperation="Compile Cna Output" fromProperty="output_dir" toOperation="Merge Cna Output By Chrom" toProperty="bam_to_cna_output_dir" />

  <link fromOperation="Merge Cna Output By Chrom" fromProperty="output_dir" toOperation="Execute" toProperty="data_directory" />
  <link fromOperation="input connector" fromProperty="data_directory" toOperation="Execute" toProperty="output_directory" />
  
  <link fromOperation="Execute" fromProperty="test_output_directory" toOperation="Individual Region Calls" toProperty="cmds_test_dir" />
  <link fromOperation="Merge Cna Output By Chrom" fromProperty="output_dir" toOperation="Individual Region Calls" toProperty="cmds_input_data_dir" />
  <link fromOperation="input connector" fromProperty="region_output_dir" toOperation="Individual Region Calls" toProperty="output_dir" />

  <link fromOperation="Individual Region Calls" fromProperty="output_dir" toOperation="Create Output Table" toProperty="region_call_dir" />
  <link fromOperation="input connector" fromProperty="table_output" toOperation="Create Output Table" toProperty="output_file" />

  <link fromOperation="Create Output Table" fromProperty="output_file" toOperation="output connector" toProperty="final_output" />
  
  <operation name="Compile Cna Output">
    <operationtype commandClass="Genome::Model::Tools::Cmds::CompileCnaOutput" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Merge Cna Output By Chrom">
    <operationtype commandClass="Genome::Model::Tools::Cmds::MergeCnaOutputByChrom" typeClass="Workflow::OperationType::Command" />
  </operation>
 

  <operation name="Execute">
    <operationtype commandClass="Genome::Model::Tools::Cmds::Execute" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Individual Region Calls">
    <operationtype commandClass="Genome::Model::Tools::Cmds::IndividualRegionCalls" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Create Output Table">
    <operationtype commandClass="Genome::Model::Tools::Cmds::CreateOutputTable" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">input_models</inputproperty>
    <inputproperty isOptional="Y">compile_cna_output_dir</inputproperty>
    <inputproperty isOptional="Y">merge_output_dir</inputproperty>
    <inputproperty isOptional="Y">region_output_dir</inputproperty>
    <inputproperty isOptional="Y">table_output</inputproperty>
    <inputproperty>data_directory</inputproperty>
    <outputproperty>final_output</outputproperty>
  </operationtype>

</workflow>


