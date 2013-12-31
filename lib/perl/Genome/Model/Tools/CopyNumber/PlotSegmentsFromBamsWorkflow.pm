package Genome::Model::Tools::CopyNumber::PlotSegmentsFromBamsWorkflow;

use warnings;
use strict;
use Genome;
use File::Basename;
use Workflow::Simple;

class Genome::Model::Tools::CopyNumber::PlotSegmentsFromBamsWorkflow {
    is => 'Command',
    has => [
        normal_bam => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc => 'Normal .bam file',
        },
        tumor_bam => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc => 'Tumor .bam file',
        },
        output_directory => {
            is => 'String',
            is_input => 1,
            doc => 'Place for results to land',
        },
    ],
    has_optional => [
        gaps_file => {
            is => 'String',
            doc => 'The gaps file to pass along to CnvSeg. If this is not provided, it will be looked up via the genome build.',
        },
        centromere_file => {
            is => 'String',
            doc => 'The centromere file to pass along to CnvSeg. If this is not provided, it will be looked up via the genome build.',
        },
        genome_build => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            valid_values => ['36', '37', 'mm9'],
            doc => "Reference sequence to use",
        },
        output_pdf => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc => 'Name of output file which is the CN graph (PDF)',
            default => 'cnv_graph.pdf',
        },
        sex => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc => "choose 'male' or 'female'",
            default => 'male'
        },
        max_copy_number => {
            is => 'String',
            is_input => 1,
            doc => 'Max copy number',
            default => 4,
        },
        plot_ymax => {
            is => 'Number',
            is_input => 1,
            is_output => 1,
            doc => 'set the max value of the y-axis on the CN plots',
            default => '6',
        },
        bam2cn_window => {
            is => 'Number',
            is_input => 1,
            is_output => 1,
            doc => 'set the window-size used for the single-genome CN estimation',
            default => '2500',
        },
        cnvseg_markers => {
            is => 'Number',
            is_input => 1,
            is_output => 1,
            doc => 'number of consecutive markers needed to make a CN gain/loss prediction',
            default => '4',
        },
        lowres => {
            is => 'String',
            is_input => 1,
            is_output => 1,
            doc => 'Set this value to zero for higher resolution output',
            default => 1,
        },
    ],
    has_param => [
        lsf_queue => {
            default => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
        }
    ],

};

sub help_detail {
    "This module takes a pair of bams and graphs their copy-number variations on a pdf for easy comparison"
}
    

sub execute {

    my $self = shift;
    my $output_dir = $self->output_directory;
    unless(-d $output_dir){
        Genome::Sys->create_directory($output_dir);
    }
    unless(-e $self->normal_bam){
        die $self->error_message("Could not locate normal bam file at: ".$self->normal_bam);
    }
    unless(-e $self->tumor_bam){
        die $self->error_message("Could not locate tumor bam file at: ".$self->tumor_bam);
    }
    unless(-d $output_dir){
        die $self->error_message("Could not locate output directory.");
    }
    my $normal_bam_name = basename($self->normal_bam);
    my $tumor_bam_name = basename($self->tumor_bam);

    $self->_resolve_gaps_and_centromere_files;

    my $workflow = Workflow::Operation->create_from_xml(\*DATA);

    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    my %input;    

    my $log_dir = $self->output_directory;
    if(Workflow::Model->parent_workflow_log_dir) {
        $log_dir = Workflow::Model->parent_workflow_log_dir;
    }

    $workflow->log_dir($log_dir);

    $input{normal_bam} =  $self->normal_bam;
    $input{normal_bamtocn} = $output_dir."/".$normal_bam_name.".bamtocn";
    $input{normal_cnv_seg} = $output_dir."/".$normal_bam_name.".cnvseg";

    $input{tumor_bam} =  $self->tumor_bam;
    $input{tumor_bamtocn} = $output_dir."/".$tumor_bam_name.".bamtocn";
    $input{tumor_cnv_seg} = $output_dir."/".$tumor_bam_name.".cnvseg";

    $input{window_size} = $self->bam2cn_window;
    $input{min_markers} = $self->cnvseg_markers;
    $input{max_copy_number} = $self->max_copy_number;
    my $pdf = $self->output_pdf;
    $input{output_pdf} = $output_dir."/".$pdf;
    $input{cnvhmm_input} = 1;
    $input{lowres} = $self->lowres;
    $input{plot_title} = "Tumor: $tumor_bam_name, Normal: $normal_bam_name";
    $input{genome_build} = $self->genome_build;
    $input{sex} = $self->sex;
    $input{y_max} = $self->plot_ymax;
    $input{centromere_file} = $self->centromere_file;
    $input{gaps_file} = $self->gaps_file;

    my $result = Workflow::Simple::run_workflow_lsf( $workflow, %input);

    # Collect and analyze results
    unless($result){
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

# Set the gaps_file and centromere_file based upon the genome_build, if needed and applicable
sub _resolve_gaps_and_centromere_files {
    my $self = shift;
    my $refseq_version = $self->genome_build;

    for my $file_type ("gaps", "centromere") {
        my $accessor = $file_type . "_file";
        unless ($self->$accessor) {
            $self->status_message("attempting to resolve dbpath for build_version '$refseq_version' and file_type '$file_type'");
            my $file = Genome::Sys->dbpath("tgi/misc-annotation/human","build$refseq_version-20130113") . "/$file_type.csv";
            unless (-e $file) {
                die $self->error_message("Could not find the $file_type file located at $file");
            }
            $self->$accessor($file);
        }
    }

    return 1;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Plot Segments From Bams">

  <link fromOperation="input connector" fromProperty="normal_bam" toOperation="normal bam-to-cn" toProperty="aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="normal_bamtocn" toOperation="normal bam-to-cn" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="window_size" toOperation="normal bam-to-cn" toProperty="window_size" />

  <link fromOperation="input connector" fromProperty="tumor_bam" toOperation="tumor bam-to-cn" toProperty="aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="tumor_bamtocn" toOperation="tumor bam-to-cn" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="window_size" toOperation="tumor bam-to-cn" toProperty="window_size" />

  <link fromOperation="normal bam-to-cn" fromProperty="output_file" toOperation="normal cnv-seg" toProperty="copy_number_file" />
  <link fromOperation="input connector" fromProperty="normal_cnv_seg" toOperation="normal cnv-seg" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="min_markers" toOperation="normal cnv-seg" toProperty="min_markers" />
  <link fromOperation="input connector" fromProperty="max_copy_number" toOperation="normal cnv-seg" toProperty="max_copy_number" />
  <link fromOperation="input connector" fromProperty="gaps_file" toOperation="normal cnv-seg" toProperty="gap_file" />
  <link fromOperation="input connector" fromProperty="centromere_file" toOperation="normal cnv-seg" toProperty="centromere_file" />

  <link fromOperation="tumor bam-to-cn" fromProperty="output_file" toOperation="tumor cnv-seg" toProperty="copy_number_file" />
  <link fromOperation="input connector" fromProperty="tumor_cnv_seg" toOperation="tumor cnv-seg" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="min_markers" toOperation="tumor cnv-seg" toProperty="min_markers" />
  <link fromOperation="input connector" fromProperty="max_copy_number" toOperation="tumor cnv-seg" toProperty="max_copy_number" />
  <link fromOperation="input connector" fromProperty="gaps_file" toOperation="tumor cnv-seg" toProperty="gap_file" />
  <link fromOperation="input connector" fromProperty="centromere_file" toOperation="tumor cnv-seg" toProperty="centromere_file" />

  <link fromOperation="tumor cnv-seg" fromProperty="output_file" toOperation="plot segments" toProperty="tumor_segment_file" />
  <link fromOperation="normal cnv-seg" fromProperty="output_file" toOperation="plot segments" toProperty="normal_segment_file" />
  <link fromOperation="input connector" fromProperty="output_pdf" toOperation="plot segments" toProperty="output_pdf" />
  <link fromOperation="input connector" fromProperty="cnvhmm_input" toOperation="plot segments" toProperty="cnvhmm_input" />
  <link fromOperation="input connector" fromProperty="lowres" toOperation="plot segments" toProperty="lowres" />
  <link fromOperation="input connector" fromProperty="plot_title" toOperation="plot segments" toProperty="plot_title" />
  <link fromOperation="input connector" fromProperty="genome_build" toOperation="plot segments" toProperty="genome_build" />
  <link fromOperation="input connector" fromProperty="sex" toOperation="plot segments" toProperty="sex" />
  <link fromOperation="input connector" fromProperty="y_max" toOperation="plot segments" toProperty="ymax" />

  <link fromOperation="plot segments" fromProperty="output_pdf" toOperation="output connector" toProperty="output" />

  <operation name="tumor bam-to-cn">
    <operationtype commandClass="Genome::Model::Tools::CopyNumber::BamToCn" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="normal bam-to-cn">
    <operationtype commandClass="Genome::Model::Tools::CopyNumber::BamToCn" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="tumor cnv-seg">
    <operationtype commandClass="Genome::Model::Tools::CopyNumber::CnvSeg" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="normal cnv-seg">
    <operationtype commandClass="Genome::Model::Tools::CopyNumber::CnvSeg" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="plot segments">
    <operationtype commandClass="Genome::Model::Tools::CopyNumber::PlotSegments" typeClass="Workflow::OperationType::Command" />
  </operation>



  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">normal_bam</inputproperty>
    <inputproperty isOptional="Y">normal_bamtocn</inputproperty>
    <inputproperty isOptional="Y">normal_cnv_seg</inputproperty>
    <inputproperty isOptional="Y">tumor_bam</inputproperty>
    <inputproperty isOptional="Y">tumor_bamtocn</inputproperty>
    <inputproperty isOptional="Y">tumor_cnv_seg</inputproperty>
    <inputproperty isOptional="Y">window_size</inputproperty>
    <inputproperty isOptional="Y">min_markers</inputproperty>
    <inputproperty isOptional="Y">max_copy_number</inputproperty>
    <inputproperty isOptional="Y">output_pdf</inputproperty>
    <inputproperty isOptional="Y">cnvhmm_input</inputproperty>
    <inputproperty isOptional="Y">lowres</inputproperty>
    <inputproperty isOptional="Y">plot_title</inputproperty>
    <inputproperty isOptional="Y">genome_build</inputproperty>
    <inputproperty isOptional="Y">sex</inputproperty>
    <inputproperty isOptional="Y">y_max</inputproperty>
    <inputproperty isOptional="Y">gaps_file</inputproperty>
    <inputproperty isOptional="Y">centromere_file</inputproperty>
    <outputproperty>output</outputproperty>
  </operationtype>

</workflow>
