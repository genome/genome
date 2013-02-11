package Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow;

use warnings;
use strict;

use Cwd;
use Genome;
use Workflow::Simple;

my $DEFAULT_VERSION = '0.1';

class Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces somatic copy-number calls from paired samples",
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>10000] span[hosts=1]',
        },
    ],
};


sub _detect_variants {
    my $self = shift;

    ##parse input params string - expected format
    #--bamwindow-version 0.4 --bamwindow-params [-w 10000 -r -l -s]
    my $params = $self->params;
    my $bamwindow_version;
    my $bamwindow_params;

    if ($params =~ m/--bamwindow-version/) {
        $bamwindow_version =~ m/--bamwindow-version\s*(\d+\.?\d?)\s*/;
        $bamwindow_params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
    }

    my $copycat_params;
    #--copycat-params [--per-read-length --per-library]
    if ($params =~ m/--copycat-params/) {
        $bamwindow_params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
    }
    my $per_library = 0;
    my $per_read_length = 0;
    if($copycat_params =~ m/--per-library/) {
        $per_library = 1;
    }
    if($copycat_params =~ m/--per-read-length/) {
        $per_read_length = 1;
    }



    my %input;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml(\*DATA);

    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    # Collect and set input parameters
    #for bam-window
    $input{bamwindow_version} = $bamwindow_version;
    $input{bamwindow_params} = $bamwindow_params;
    $input{tumor_bam_file} = $self->aligned_reads_input;
    $input{normal_bam_file} = $self->control_aligned_reads_input;
    #for copycat
    $input{copycat_params} = $copycat_params;
    $input{per_read_length} = $per_read_length;
    $input{per_library} = $per_library;
    $input{tumor_samtools_file} = get_samtools_results();
    $input{reference_build_id} = $self->reference_build_id;


    $self->_dump_workflow($workflow);

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    # Launch workflow
    $self->status_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input);

    # Collect and analyze results
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }
    $self->_workflow_result($result);

    return 1;
}



sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

sub _sort_detector_output {
    return 0;
}


##todo comment
sub get_samtools_results{
    my ($self, $alignment_result) = @_;
    my $return_snp_file;
    my $bam_path = $self->aligned_reads_input;
    my @results = Genome::Model::Tools::DetectVariants2::Result::Detector->get(
        detector_name => "Genome::Model::Tools::DetectVariants2::Samtools",
        aligned_reads =>$bam_path );
    if(@results) {
        return $results[0]->path("snvs.hq");
    } else {
        $self->status_message("Could not find any DV2::Samtools result object for $bam_path");
        #alternative lookup - maybe later?
    }
    return "NULL";
}


1;

__DATA__
<?xml version='1.0' standalone='yes'?>

    <workflow name="CopyCatSomatic Detect Variants Module">


    <link fromOperation="input connector" fromProperty="normal_bam_file" toOperation="BamWindow Normal" toProperty="aligned_reads_input"/>
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Normal" toProperty="params" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Normal" toProperty="version" />


    <link fromOperation="input connector" fromProperty="tumor_bam_file" toOperation="BamWindow Tumor" toProperty="aligned_reads_input" />
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Tumor" toProperty="params" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Tumor" toProperty="version" />


    <link fromOperation="BamWindow Normal" fromProperty="output_directory" toOperation="CopyCatSomatic" toProperty="normal_window_dir" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_directory" toOperation="CopyCatSomatic" toProperty="tumor_window_dir" />
    <link fromOperation="input connector" fromProperty="per_library" toOperation="CopyCat Somatic" toProperty="per_library" />
    <link fromOperation="input connector" fromProperty="per_readlength" toOperation="CopyCat Somatic" toProperty="per_readlength" />
    <link fromOperation="input connector" fromProperty="tumor_samtools_file" toOperation="CopyCat Somatic" toProperty="tumor_samtools_file" />
    <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="CopyCat Somatic" toProperty="reference_build_id" />
    <link fromOperation="input connector" fromProperty="copycat_params" toOperation="CopyCat Somatic" toProperty="params" />


    <link fromOperation="BamWindow Normal" fromProperty="output_directory" toOperation="output connector" toProperty="bam_window_normal_output_directory" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_directory" toOperation="output connector" toProperty="bam_window_tumor_output_directory" />
    <link fromOperation="CopyCat Somatic" fromProperty="output_directory" toOperation="output connector" toProperty="copycat_output_directory" />

  <operation name="BamWindow Normal">
    <operationtype commandClass="Genome::Model::Tools::DetectVariants2::BamWindow" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="BamWindow Tumor">
    <operationtype commandClass="Genome::Model::Tools::DetectVariants2::BamWindow" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="CopyCat Somatic">
    <operationtype commandClass="Genome::Model::Tools::DetectVariants2::CopyCatSomatic" typeClass="Workflow::OperationType::Command" />
  </operation>


  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>tumor_bam_file</inputproperty>
    <inputproperty>normal_bam_file</inputproperty>
    <inputproperty>bamwindow_params</inputproperty>
    <inputproperty>bamwindow_version</inputproperty>
    <inputproperty>reference_build_id</inputproperty>
    <inputproperty>per_read_length</inputproperty>
    <inputproperty>per_library</inputproperty>
    <inputproperty>tumor_samtools_file</inputproperty>

    <outputproperty>bam_window_normal_output_directory</outputproperty>
    <outputproperty>bam_window_tumor_output_directory</outputproperty>
    <outputproperty>copycat_output_directory</outputproperty>
  </operationtype>

</workflow>
