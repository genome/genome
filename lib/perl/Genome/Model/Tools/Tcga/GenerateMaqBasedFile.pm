package Genome::Model::Tools::Tcga::GenerateMaqBasedFile;

use strict;
use warnings;

use Genome;
use Workflow;
 
class Genome::Model::Tools::Tcga::GenerateMaqBasedFile {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); },
    has => [
        workflow_log_directory => {
                    is => 'String',
                    doc => 'The directory where the workflow logs (LSF output) should be dumped.' ,
        },
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

    my $model_id = $self->model_id;

    $self->_operation->log_dir($self->workflow_log_directory);

    $self->status_message("Generating TCGA compliant file for model id:".$model_id);
    $self->status_message("Using working directory:".$working_dir);
    $self->status_message("Workflow log directory:".$self->workflow_log_directory);
    $self->status_message("Delete intermediate files on completion: ".$self->cleanup);
    $self->status_message("Creating required directories.");

    Genome::Sys->create_directory("$working_dir");
    Genome::Sys->create_directory("$working_dir/header");
    Genome::Sys->create_directory("$working_dir/aligned");
    Genome::Sys->create_directory("$working_dir/unaligned");
    Genome::Sys->create_directory("$working_dir/maps");
    Genome::Sys->create_directory("$working_dir/logs");
    Genome::Sys->create_directory("$working_dir/tmp");
    Genome::Sys->create_directory("$working_dir/sams");
    Genome::Sys->create_directory("$working_dir/bams");

    $self->status_message("Pre-execute of GenerateMaqBasedFile complete.");

    return 1;
}

sub post_execute {
    my $self = shift;
    my $working_dir = $self->working_directory;

    my $cleanup = $self->cleanup;

    if ($cleanup) {
        $self->status_message("Cleaning up intermediate files.");
        unlink(<$working_dir/header/*>);
        unlink(<$working_dir/aligned/*>);
        unlink(<$working_dir/unaligned/*>);
        unlink(<$working_dir/maps/*>);
        unlink(<$working_dir/tmp/*>);
        unlink(<$working_dir/sams/*>);
        unlink(<$working_dir/bams/*>);
    } else {
        $self->status_message("Leaving intermediate files behind. Done.");
    }

    return 1;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Generate Maq Based Tcga File">

  <link fromOperation="input connector" fromProperty="model_id"                     toOperation="Maq Unaligned" toProperty="model_id" />
  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="Maq Unaligned" toProperty="working_directory" />
  <link fromOperation="Maq Unaligned"   fromProperty="unaligned_sam_file_directory" toOperation="Merge All"     toProperty="unaligned_sam_file_directory" />
  <link fromOperation="Maq Unaligned"   fromProperty="read_group_directory"         toOperation="Merge All"     toProperty="read_group_directory" />
  <link fromOperation="Maq Unaligned"   fromProperty="program_group_directory"      toOperation="Merge All"     toProperty="program_group_directory" />
  
  <link fromOperation="input connector" fromProperty="model_id"                     toOperation="Maq Aligned"   toProperty="model_id" />
  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="Maq Aligned"   toProperty="working_directory" />
  <link fromOperation="Maq Aligned"     fromProperty="aligned_sam_file_directory"   toOperation="Merge All"     toProperty="aligned_sam_file_directory" />
  
  <link fromOperation="input connector" fromProperty="model_id"                     toOperation="Merge All"     toProperty="model_id" />
  <link fromOperation="input connector" fromProperty="working_directory"            toOperation="Merge All"     toProperty="working_directory" />
  <link fromOperation="Merge All"       fromProperty="final_file"                   toOperation="output connector" toProperty="final_file" />


<operation name="Maq Unaligned">
    <operationtype commandClass="Genome::Model::Tools::Tcga::ConvertUnalignedMaqReadsToSamFiles" typeClass="Workflow::OperationType::Command" />
</operation>
<operation name="Maq Aligned">
    <operationtype commandClass="Genome::Model::Tools::Tcga::ConvertAlignedMapsToSamFiles" typeClass="Workflow::OperationType::Command" />
</operation>
<operation name="Merge All">
    <operationtype commandClass="Genome::Model::Tools::Tcga::MergePerLaneSamFiles" typeClass="Workflow::OperationType::Command" />
</operation>

<operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>model_id</inputproperty>
    <inputproperty>working_directory</inputproperty>
    <outputproperty>final_file</outputproperty>
</operationtype>

</workflow>
