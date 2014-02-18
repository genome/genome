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
    has => [
        bamwindow_filter_to_chromosomes => {
            is => 'Text',
            is_optional => 1,
            is_many => 1,
            doc => 'Chromosomes that bamwindow will filter output to (generally for test purposes)',
        },
        annotation_version => {
            is => 'Text',
            doc => 'version of the copycat annotation to use',
            is_optional => 1,
        }
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
        $params =~ m/--bamwindow-version\s*(\d+\.?\d?)\s*/;
        $bamwindow_version = $1;
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params  = $1;
        $bamwindow_params =~ s/--bamwindow-params\s*//;
    }

    my $copycat_params;
    #--copycat-params [--per-read-length --per-library]
    if ($params =~ m/--copycat-params/) {
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params = $1;
        $params =~ m/--copycat-params\s*\{([^\}]+)\}/;
        $copycat_params = $1;
    }
    my $per_library = 0;
    my $per_read_length = 0;
    if($copycat_params =~ m/--per-library/) {
        $per_library = 1;
    }
    if($copycat_params =~ m/--per-read-length/) {
        $per_read_length = 1;
    }

    #annotation-version
    my $annotation_version;
    if($params =~ m/--annotation-version/){
        $params =~ m/--annotation-version\s*(\d+\.?\d?)\s*/;
        $annotation_version = $1;
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
    $input{tumor_window_file} = join('/', $self->_temp_staging_directory, 'tumor_window_file');
    $input{normal_window_file} = join('/', $self->_temp_staging_directory, 'normal_window_file');
    my @bamwindow_filter_to_chromosomes = $self->bamwindow_filter_to_chromosomes;
    $input{bamwindow_filter_to_chromosomes} = \@bamwindow_filter_to_chromosomes;

    #for copycat
    $input{per_read_length} = $per_read_length;
    $input{per_library} = $per_library;
    $input{tumor_samtools_file} = $self->get_samtools_results($self->aligned_reads_input);
    $input{normal_samtools_file} = $self->get_samtools_results($self->control_aligned_reads_input);
    $input{copycat_output_directory} = $self->_temp_staging_directory;
    $input{annotation_version} = $annotation_version;
    $input{reference_build_id} = $self->reference_build_id;

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }
    $workflow->log_dir($log_dir);

    # Launch workflow
    $self->debug_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %input);

    # Collect and analyze results
    unless($result){
        if (@Workflow::Simple::ERROR){
            print Data::Dumper->Dumper(@Workflow::Simple::ERROR), "\n";
        }
        die $self->error_message("Workflow did not return correctly");
    }

    return 1;
}

sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

sub _sort_detector_output {
    return 1;
}

sub get_samtools_results{
    my $self = shift;
    my $bam_path = shift;
    my $return_snp_file;
    my @results = Genome::Model::Tools::DetectVariants2::Result::DetectionBase->get(
        detector_name => "Genome::Model::Tools::DetectVariants2::Samtools",
        aligned_reads =>$bam_path );
    if(@results) {
        return $results[0]->path("snvs.hq");
    } else {
        $self->debug_message("Could not find any DV2::Samtools result object for $bam_path");
        #alternative lookup - maybe later?
    }
    return "";
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>
    <workflow name="CopyCatSomatic Detect Variants Module">


    <link fromOperation="input connector" fromProperty="normal_bam_file" toOperation="BamWindow Normal" toProperty="bam_file"/>
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Normal" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Normal" toProperty="version" />
    <link fromOperation="input connector" fromProperty="normal_window_file" toOperation="BamWindow Normal" toProperty="output_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_filter_to_chromosomes" toOperation="BamWindow Normal" toProperty="filter_to_chromosomes" />


    <link fromOperation="input connector" fromProperty="tumor_bam_file" toOperation="BamWindow Tumor" toProperty="bam_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Tumor" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Tumor" toProperty="version" />
    <link fromOperation="input connector" fromProperty="tumor_window_file" toOperation="BamWindow Tumor" toProperty="output_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_filter_to_chromosomes" toOperation="BamWindow Tumor" toProperty="filter_to_chromosomes" />


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="normal_window_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="tumor_window_file" />
    <link fromOperation="input connector" fromProperty="per_library" toOperation="CopyCat Somatic" toProperty="per_library" />
    <link fromOperation="input connector" fromProperty="per_read_length" toOperation="CopyCat Somatic" toProperty="per_read_length" />
    <link fromOperation="input connector" fromProperty="tumor_samtools_file" toOperation="CopyCat Somatic" toProperty="tumor_samtools_file" />
    <link fromOperation="input connector" fromProperty="normal_samtools_file" toOperation="CopyCat Somatic" toProperty="normal_samtools_file" />
    <link fromOperation="input connector" fromProperty="copycat_output_directory" toOperation="CopyCat Somatic" toProperty="output_directory" />
    <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="CopyCat Somatic" toProperty="reference_build_id" />
    <link fromOperation="input connector" fromProperty="annotation_version" toOperation="CopyCat Somatic" toProperty="annotation_version" />


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_normal_output_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_tumor_output_file" />
    <link fromOperation="CopyCat Somatic" fromProperty="output_directory" toOperation="output connector" toProperty="copycat_output_directory" />

    <operation name="BamWindow Normal">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="BamWindow Tumor">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="CopyCat Somatic">
    <operationtype commandClass="Genome::Model::Tools::CopyCat::Somatic" typeClass="Workflow::OperationType::Command" />
    </operation>


    <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>tumor_bam_file</inputproperty>
    <inputproperty>normal_bam_file</inputproperty>
    <inputproperty>tumor_window_file</inputproperty>
    <inputproperty>normal_window_file</inputproperty>
    <inputproperty>bamwindow_params</inputproperty>
    <inputproperty>bamwindow_version</inputproperty>
    <inputproperty>per_read_length</inputproperty>
    <inputproperty>per_library</inputproperty>
    <inputproperty>tumor_samtools_file</inputproperty>
    <inputproperty>normal_samtools_file</inputproperty>
    <inputproperty>copycat_output_directory</inputproperty>
    <inputproperty>annotation_version</inputproperty>
    <inputproperty>reference_build_id</inputproperty>
    <inputproperty isOptional="Y">bamwindow_filter_to_chromosomes</inputproperty>

    <outputproperty>bam_window_normal_output_file</outputproperty>
    <outputproperty>bam_window_tumor_output_file</outputproperty>
    <outputproperty>copycat_output_directory</outputproperty>
    </operationtype>

    </workflow>
