package Genome::Model::Somatic::Command::RunPindel;

use strict;
use warnings;

use Workflow;
use File::Basename;

class Genome::Model::Somatic::Command::RunPindel {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); }
};

sub help_brief {
    "Runs the pindel pipeline on the last complete build of a somatic model."
}

sub help_synopsis{
    my $self = shift;
    return <<"EOS"
genome model somatic pindel --model-id 123 --output-directory /someplace/for/output (do not put this in an allocated build directory, it will make allocations inaccurate)
EOS
}

sub help_detail {
    my $self = shift;
    return <<"EOS"
This tool runs the pindel pipeline on the last complete build of a somatic model. This pipeline will be integrated with the somatic pipeline in the future.
EOS
}

sub pre_execute {
    my $self = shift;
    # Obtain normal and tumor bams and check them. Either from somatic model id or from direct specification. 
    my ($build, $tumor_bam, $normal_bam);
    if ( ($self->model_id) && ($self->tumor_bam || $self->normal_bam) ) {
        $self->error_message("Usage error. Please specify either model_id OR tumor_bam and normal_bam, not both");
        die;
    } elsif ($self->model_id) {
        my $model = Genome::Model::Somatic->get($self->model_id);
        unless ($model) {
            $self->status_message("Could not locate somatic model for model-id ".$self->model_id." checking somatic-capture models.");
            $model = Genome::Model::SomaticCapture->get($self->model_id);
            unless($model){
                $self->error_message("Could not get a somatic capture model for id " . $self->model_id);
                die;
            }
        }
        my $tumor_ra_model = $model->tumor_model;
        my $normal_ra_model = $model->normal_model;
        unless ($tumor_ra_model && $normal_ra_model) {
            $self->error_message("Could not get tumor or normal model from somatic model " . $self->model_id);
            die;
        }
        my $normal_ra_build = $normal_ra_model->last_succeeded_build;
        my $tumor_ra_build = $tumor_ra_model->last_succeeded_build;
        unless ($tumor_ra_build && $normal_ra_build) {
            $self->error_message("Could not get tumor or normal model from somatic model " . $self->model_id);
            die;
        }

        $normal_bam = $normal_ra_build->whole_rmdup_bam_file;
        $tumor_bam = $tumor_ra_build->whole_rmdup_bam_file;
        unless(-s $tumor_bam && -s $normal_bam){
            $self->error_message("Couldn't locate tumor or normal bam files.");
            die;
        }
        #stuff the bam paths back into self so they are accessible to the workflow input connector.
        $self->tumor_bam($tumor_bam);
        $self->normal_bam($normal_bam);
    } elsif ($self->tumor_bam || $self->normal_bam) {
        if(defined($self->normal_bam)){
            $normal_bam = $self->normal_bam;
        }
        $tumor_bam = $self->tumor_bam;
    } else {
        $self->error_message("Usage error. Please specify either model_id OR tumor_bam and normal_bam");
        die;
    }
    
    unless (-s $normal_bam || not defined($normal_bam)) {
        $self->error_message("Normal bam $normal_bam does not exist or has 0 size");
        die;
    }

    unless (-s $tumor_bam) {
        $self->error_message("tumor bam $tumor_bam does not exist or has 0 size");
        die;
    }
    unless(defined($self->reference_sequence_input)){
        $self->reference_sequence_input( Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa' );
    }
    unless(defined($self->tier_file_directory)){
        $self->tier_file_directory('/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/');
    }

    # Set default params
    unless ($self->annotate_no_headers) { $self->annotate_no_headers(1); }
    unless ($self->transcript_annotation_filter) { $self->transcript_annotation_filter("top"); }

    if($self->chromosome_list) { 
        my @chrom_list = split ",", $self->chromosome_list;
        $self->chromosome_list(\@chrom_list);
    }
    else {
        my @chrom_list = $self->get_chromosome_list($self->reference_sequence_input);
        $self->chromosome_list(\@chrom_list); 
    }
    unless ($self->indel_bed_output) { $self->indel_bed_output($self->output_directory . '/indels_all_sequences.bed'); }

    unless(defined($self->version)){
        $self->version('0.2');
    }
    unless($self->use_old_pindel){
        $self->use_old_pindel(1);
    }

    my %default_filenames = $self->default_filenames;
    for my $param (keys %default_filenames) {
        # set a default param if one has not been specified
        my $default_filename = $default_filenames{$param};
        $self->$param( join('/', $self->output_directory, $default_filename) );
    }

    # create directories
    #for my $directory ( $self->assemble_t1n_dir, $self->assemble_t1t_dir, $self->assemble_t2n_dir, $self->assemble_t2t_dir, $self->assemble_t3n_dir, $self->assemble_t3t_dir) {
    #    unless ( Genome::Sys->create_directory($directory) ) {
    #        $self->error_message("Failed to create directory $directory");
    #        die;
    #    }
    #}

    return 1;
}

sub get_chromosome_list {
    my $self = shift;
    my $fasta = shift;
    my ($filename, $path,$suffix) = fileparse($fasta,['.fa','.fasta']);
    my $index_file = $path."/all_sequences.fasta.fai";
    unless(-e $index_file){
        $self->error_message("Could not locate reference sequence index file at ".$index_file);
        die $self->error_message;
    }
    my $fh = Genome::Sys->open_file_for_reading($index_file);
    my @chrom_list;
    while ( my $line = $fh->getline){
        chomp $line;
        my ($chr) = split "\t", $line;

        #for now, we are not including MT or NT contigs
        if($chr =~ m/[MT,NT]/i){
            next;
        }
        push @chrom_list, $chr;
    }

    return @chrom_list;
}

sub default_filenames{
    my $self = shift;
   
    my %default_filenames = (
        annotation_output => "tier1_annotated.csv",
        intersect_output => "somatic_intersected.bed",
    );

    return %default_filenames;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Pindel No Assembly" logDir="/gscmnt/ams1158/info/pindel/logs/pindel_no_assembly">

  <link fromOperation="input connector" fromProperty="normal_bam" toOperation="Pindel" toProperty="control_aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="tumor_bam" toOperation="Pindel" toProperty="aligned_reads_input" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Pindel" toProperty="output_directory" />
  <link fromOperation="input connector" fromProperty="chromosome_list" toOperation="Pindel" toProperty="chromosome" />
  <link fromOperation="input connector" fromProperty="version" toOperation="Pindel" toProperty="version" />
  <link fromOperation="input connector" fromProperty="reference_sequence_input" toOperation="Pindel" toProperty="reference_sequence_input" />

  <link fromOperation="Pindel" fromProperty="indel_bed_output" toOperation="Cat" toProperty="source" />
  <link fromOperation="input connector" fromProperty="indel_bed_output" toOperation="Cat" toProperty="dest" />

  <link fromOperation="Cat" fromProperty="dest" toOperation="Pre-Assembly Tiering" toProperty="variant_bed_file" />
  <link fromOperation="input connector" fromProperty="tier_file_directory" toOperation="Pre-Assembly Tiering" toProperty="tier_file_location" />

  <link fromOperation="Pre-Assembly Tiering" fromProperty="tier1_output" toOperation="Annotation" toProperty="variant_bed_file" />
  <link fromOperation="input connector" fromProperty="annotation_output" toOperation="Annotation" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="annotate_no_headers" toOperation="Annotation" toProperty="no_headers" />
  <link fromOperation="input connector" fromProperty="transcript_annotation_filter" toOperation="Annotation" toProperty="annotation_filter" />
  

  <link fromOperation="Pre-Assembly Tiering" fromProperty="tier1_output" toOperation="Calculate Read Support Tier1" toProperty="indels_all_sequences_bed_file" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Calculate Read Support Tier1" toProperty="pindel_output_directory" />
  <link fromOperation="input connector" fromProperty="use_old_pindel" toOperation="Calculate Read Support Tier1" toProperty="use_old_pindel" />
  <link fromOperation="Pre-Assembly Tiering" fromProperty="tier2_output" toOperation="Calculate Read Support Tier2" toProperty="indels_all_sequences_bed_file" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Calculate Read Support Tier2" toProperty="pindel_output_directory" />
  <link fromOperation="input connector" fromProperty="use_old_pindel" toOperation="Calculate Read Support Tier2" toProperty="use_old_pindel" />
  <link fromOperation="Pre-Assembly Tiering" fromProperty="tier3_output" toOperation="Calculate Read Support Tier3" toProperty="indels_all_sequences_bed_file" />
  <link fromOperation="input connector" fromProperty="output_directory" toOperation="Calculate Read Support Tier3" toProperty="pindel_output_directory" />
  <link fromOperation="input connector" fromProperty="use_old_pindel" toOperation="Calculate Read Support Tier3" toProperty="use_old_pindel" />

  <link fromOperation="Calculate Read Support Tier1" fromProperty="_output_filename" toOperation="Filter Read Support Tier1" toProperty="read_support_file" />
  <link fromOperation="Filter Read Support Tier1" fromProperty="output_file" toOperation="output connector" toProperty="tier_1_filtered_read_support" />
  <link fromOperation="Calculate Read Support Tier2" fromProperty="_output_filename" toOperation="Filter Read Support Tier2" toProperty="read_support_file" />
  <link fromOperation="Filter Read Support Tier2" fromProperty="output_file" toOperation="output connector" toProperty="tier_2_filtered_read_support" />
  <link fromOperation="Calculate Read Support Tier3" fromProperty="_output_filename" toOperation="Filter Read Support Tier3" toProperty="read_support_file" />
  <link fromOperation="Filter Read Support Tier3" fromProperty="output_file" toOperation="output connector" toProperty="tier_3_filtered_read_support" />

  <link fromOperation="Annotation" fromProperty="output_file" toOperation="output connector" toProperty="output" />
  
  <operation name="Pindel" parallelBy="chromosome">
    <operationtype commandClass="Genome::Model::Tools::DetectVariants::Somatic::Pindel" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Cat">
    <operationtype commandClass="Genome::Model::Tools::Cat" typeClass = "Workflow::OperationType::Command" />
  </operation>

  <operation name="Pre-Assembly Tiering">
    <operationtype commandClass="Genome::Model::Tools::FastTier::FastTier" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Calculate Read Support Tier1">
    <operationtype commandClass="Genome::Model::Tools::Somatic::CalculatePindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Calculate Read Support Tier2">
    <operationtype commandClass="Genome::Model::Tools::Somatic::CalculatePindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Calculate Read Support Tier3">
    <operationtype commandClass="Genome::Model::Tools::Somatic::CalculatePindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Filter Read Support Tier1">
    <operationtype commandClass="Genome::Model::Tools::Somatic::FilterPindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Filter Read Support Tier2">
    <operationtype commandClass="Genome::Model::Tools::Somatic::FilterPindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Filter Read Support Tier3">
    <operationtype commandClass="Genome::Model::Tools::Somatic::FilterPindelReadSupport" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Annotation">
    <operationtype commandClass="Genome::Model::Tools::Annotate::TranscriptVariants" typeClass="Workflow::OperationType::Command" />
  </operation>
  
  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty isOptional="Y">model_id</inputproperty>
    <inputproperty isOptional="Y">normal_bam</inputproperty>
    <inputproperty isOptional="Y">tumor_bam</inputproperty>
    <inputproperty isOptional="Y">output_directory</inputproperty>
    <inputproperty isOptional="Y">version</inputproperty>
    <inputproperty isOptional="Y">annotate_no_headers</inputproperty>
    <inputproperty isOptional="Y">transcript_annotation_filter</inputproperty>
    <inputproperty isOptional="Y">annotation_output</inputproperty>
    <inputproperty isOptional="Y">intersect_output</inputproperty>
    <inputproperty isOptional="Y">chromosome_list</inputproperty>
    <inputproperty isOptional="Y">indel_bed_output</inputproperty>
    <inputproperty isOptional="Y">tiered_bed_files</inputproperty>
    <inputproperty isOptional="Y">use_old_pindel</inputproperty>
    <inputproperty isOptional="Y">reference_sequence_input</inputproperty>
    <inputproperty isOptional="Y">tier_file_directory</inputproperty>

    <outputproperty>output</outputproperty>
    <outputproperty>tier_1_filtered_read_support</outputproperty>
    <outputproperty>tier_2_filtered_read_support</outputproperty>
    <outputproperty>tier_3_filtered_read_support</outputproperty>
    
  </operationtype>

</workflow>
