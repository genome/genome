package Genome::Model::Tools::DetectVariants2::MethRatio;

use warnings;
use strict;

use Genome;
use Workflow;
use File::Copy;
use Workflow::Simple;
use Cwd;

my $DEFAULT_VERSION = '2.6';

class Genome::Model::Tools::DetectVariants2::MethRatio {
    is => ['Genome::Model::Tools::DetectVariants2::WorkflowDetectorBase'],
    doc => "Runs methyl counting script on a bsmap alignment model.",
    has => [
        chromosome_list => {
            is => 'ARRAY',
            is_optional => 1,
            doc => 'list of chromosomes to run on. This will run on a default set of chromosomes from the reference sequence if not set.',
        },
   ],
    has_transient_optional => [
        _workflow_result => {
            doc => 'Result of the workflow',
        },
        _snv_output_dir => {
            is => 'String',
            doc => 'The location of the snvs.hq file',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => 'workflow'
        },
    ],
};

sub _detect_variants {
    my $self = shift;

    # Set the chromosome_list if it is not already
    unless ($self->chromosome_list) {
        my @default_chromosomes = $self->default_chromosomes;
        $self->chromosome_list(\@default_chromosomes);
    }

    # Obtain bam and check it.
    my ($build, $bam_file);
    $bam_file = $self->aligned_reads_input;

    # Set default params
    unless ($self->snv_output) { 
        $self->snv_output($self->_temp_staging_directory. '/snvs.hq');  #???
    }

    #get ref seq fasta
    my $refbuild_id = $self->reference_build_id;
    unless($refbuild_id){
        die $self->error_message("Received no reference build id.");
    }
    print "refbuild_id = ".$refbuild_id."\n";
    my $ref_seq_build = Genome::Model::Build->get($refbuild_id);
    my $reference_fasta = $ref_seq_build->full_consensus_path('fa');


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
    $input{chromosome_list} = $self->chromosome_list;
    $input{reference_build_fasta} = $reference_fasta;
    $input{bam_file} = $self->aligned_reads_input;
    $input{output_directory}  =  $self->output_directory;
    $input{version}  =  $self->version;
    

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

sub _create_temp_directories {
    my $self = shift;
    $self->_temp_staging_directory($self->output_directory);
    $self->_temp_scratch_directory($self->output_directory);
    return 1;

    return $self->SUPER::_create_temp_directories(@_);
}

sub _generate_standard_files {
    my $self = shift;
    my $staging_dir = $self->_temp_staging_directory;
    my $output_dir  = $self->output_directory;
    my $chrom_list = $self->chromosome_list;
    my $raw_output_file = $output_dir."/snvs.hq";
    my @raw_inputs = map { $output_dir."/".$_."/snvs.hq" } @$chrom_list;
    my $cat_raw = Genome::Model::Tools::Cat->create( dest => $raw_output_file, source => \@raw_inputs);
    unless($cat_raw->execute){
        die $self->error_message("Cat command failed to execute.");
    }
    $self->SUPER::_generate_standard_files(@_);
    return 1;
}

sub _promote_staged_data {
    my $self = shift;
    return 1;
}

sub _sort_detector_output {
    my $self = shift;
    return 1;
}

sub has_version {
    my $self = shift;
    my $version = shift;
    unless(defined($version)){
        $version = $self->version;
    }
    my @versions = Genome::Model::Tools::Bsmap::MethRatioWorkflow->available_methratio_versions;

    for my $v (@versions){
        if($v eq $version){
            return 1;
        }
    }

    return 0;
}

sub default_chromosomes {
    my $self = shift;
    my $refbuild = Genome::Model::Build::ReferenceSequence->get($self->reference_build_id);
    die unless $refbuild;
    my $chromosome_array_ref = $refbuild->chromosome_array_ref;
    return $self->sort_chromosomes($refbuild->chromosome_array_ref);
}

sub default_chromosomes_as_string {
    return join(',', $_[0]->default_chromosomes);
}

sub chromosome_list_as_string {
    my $self = shift;
    my $chromosome_list = $self->chromosome_list;
    return join(',', @$chromosome_list);
}

sub params_for_detector_result {
    my $self = shift;
    my ($params) = $self->SUPER::params_for_detector_result;

    if ($self->chromosome_list) {
        $params->{chromosome_list} = $self->chromosome_list_as_string;
    } else {
        $params->{chromosome_list} = $self->default_chromosomes_as_string;
    }

    return $params;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="MethRatio Detect Variants Module">

  <link fromOperation="input connector" fromProperty="bam_file" toOperation="MethRatio" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="MethRatio" toProperty="output_directory" />
  <link fromOperation="input connector" fromProperty="chromosome_list" toOperation="MethRatio" toProperty="chromosome" />
  <link fromOperation="input connector" fromProperty="reference_build_fasta" toOperation="MethRatio" toProperty="reference" />
  <link fromOperation="input connector" fromProperty="version" toOperation="MethRatio" toProperty="version" />


  <link fromOperation="MethRatio" fromProperty="output_directory" toOperation="output connector" toProperty="output" />

  <operation name="MethRatio" parallelBy="chromosome">
    <operationtype commandClass="Genome::Model::Tools::Bsmap::MethRatioWorkflow" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">bam_file</inputproperty>
    <inputproperty isOptional="Y">output_directory</inputproperty>
    <inputproperty isOptional="Y">chromosome_list</inputproperty>
    <inputproperty isOptional="Y">version</inputproperty>
    <inputproperty isOptional="Y">reference_build_fasta</inputproperty>

    <outputproperty>output</outputproperty>
  </operationtype>

</workflow>
