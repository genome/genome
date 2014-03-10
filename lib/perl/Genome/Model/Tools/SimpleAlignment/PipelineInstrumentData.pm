package Genome::Model::Tools::SimpleAlignment::PipelineInstrumentData;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::SimpleAlignment::PipelineInstrumentData {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); },
    has => [
        cleanup => { 
                    is => 'Boolean',
                    is_optional => '1',
                    default_value => '0',
                    doc => 'A clean up flag.  Will remove intermediate files if set. Default = 0, no cleanup.',
        },
    ]
};

sub help_synopsis{
    my $self = shift;
    return "TBD";
}

sub pre_execute {
    my $self = shift;

    #make required directories if they don't exist
    my $working_dir = $self->working_directory;

    $self->_operation->log_dir($self->working_directory);

    $self->debug_message("Launching Pipeline with InstrumentData.");
    $self->debug_message("Using working directory:".$working_dir);
    $self->debug_message("Workflow log directory:".$self->working_directory);
    $self->debug_message("Delete intermediate files on completion: ".$self->cleanup);
    
    $self->debug_message("Instrument data ids: ".$self->instrument_data_id);
    my @id_files = split(/,/ , $self->instrument_data_id);
    my $list_string = join("\n",@id_files);
    $self->debug_message("Listing instrument data ids: \n".$list_string); 
    #$self->debug_message("Creating required directories.");
    $self->instrument_data_id(\@id_files);

    $self->debug_message("Pre-execute of Pipeline complete.");

    return 1;
}

sub post_execute {
    my $self = shift;
    my $working_dir = $self->working_directory;

    my $cleanup = $self->cleanup;

    if ($cleanup) {
        $self->debug_message("Cleaning up intermediate files.");
     
    } else {
        $self->debug_message("Leaving intermediate files behind.");
    }

	$self->debug_message("Done.");
    return 1;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="HMP Metagenomic Pipeline for Instrument Data">

  <link fromOperation="input connector" fromProperty="working_directory"	toOperation="Align" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="instrument_data_id"	toOperation="Align" toProperty="instrument_data_id" /> 
  <link fromOperation="input connector" fromProperty="reference_name"	        toOperation="Align" toProperty="reference_name" /> 
  <link fromOperation="Align" fromProperty="aligned_file"	    		toOperation="MergeAlignments" toProperty="alignment_files" /> 

  <link fromOperation="input connector" fromProperty="working_directory"    toOperation="MergeAlignments" toProperty="working_directory" /> 
  <link fromOperation="input connector" fromProperty="merged_file"          toOperation="MergeAlignments" toProperty="merged_file" /> 
  <link fromOperation="MergeAlignments"	fromProperty="merged_file"          toOperation="output connector" toProperty="merged_file_path" />
  
  
<operation name="Align" parallelBy="instrument_data_id">
    <operationtype commandClass="Genome::Model::Tools::SimpleAlignment::AlignWrapper" typeClass="Workflow::OperationType::Command">
    </operationtype>
</operation>

<operation name="MergeAlignments">
    <operationtype commandClass="Genome::Model::Tools::SimpleAlignment::MergeAlignments" typeClass="Workflow::OperationType::Command" />
</operation>


<operationtype typeClass="Workflow::OperationType::Model">
	<inputproperty>working_directory</inputproperty>	
	<inputproperty>merged_file</inputproperty>
	<inputproperty>instrument_data_id</inputproperty>
	<inputproperty>reference_name</inputproperty>
    <outputproperty>merged_file_path</outputproperty>
</operationtype>

</workflow>
