package Genome::Model::Tools::Germline::CaptureBams;

########################################################################################################################
# CaptureBams.pm - A module for comparing tumor-normal BAM files in capture data
#
#   
########################################################################################################################

use strict;
use warnings;

use Workflow;

class Genome::Model::Tools::Germline::CaptureBams {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); }
};

sub help_brief {
    "Runs the **capture** germline pipeline workflow."
}

sub help_synopsis{
    my $self = shift;
    return <<"EOS"

example:
gmt germline capture-bams --build-id=101625141 --filtered-indelpe-snps='/gscmnt/sata835/info/medseq/model_data/2852971605/build101625141/sam_snp_related_metrics/filtered.indelpe.snps' --indels-all-sequences-filtered='/gscmnt/sata835/info/medseq/model_data/2852971605/build101625141/sam_snp_related_metrics/indels_all_sequences.filtered' --germline-bam-file='/gscmnt/sata835/info/medseq/model_data/2852971605/build101625141/alignments/101625141_merged_rmdup.bam' --data-directory=/gscmnt/sata424/info/medseq/Freimer-Boehnke/Germline_Pipeline_Test/ --regions-file=/gscmnt/sata424/info/medseq/Freimer-Boehnke/targets/capture-targets.set1.tsv

EOS
}

sub help_detail {
    my $self = shift;
    return <<"EOS"
This tool runs the capture germline pipeline to take in samtools SNPs and indels and a bam file. It results in running varscan, annotating, then outputting tiered SNPs and indels.
EOS
}

sub pre_execute {
    my $self = shift;

    # If data directory was provided... make sure it exists and set all of the file names
    if ($self->data_directory) {
        unless (-d $self->data_directory) {
            $self->error_message("Data directory " . $self->data_directory . " does not exist. Please create it.");
            return 0;
        }
        
        my %default_filenames = $self->default_filenames;
        for my $param (keys %default_filenames) {
            # set a default param if one has not been specified
            my $default_filename = $default_filenames{$param};
            unless ($self->$param) {
                #$self->status_message("Param $param was not provided... generated $default_filename as a default");
                $self->$param($self->data_directory . "/$default_filename");
            }
        }
    }

    unless (defined $self->skip_if_output_present) {
        $self->skip_if_output_present(1);
    }

    # Dummy regions_file
    unless (defined $self->regions_file) {
        $self->regions_file('/gscmnt/sata424/info/medseq/Freimer-Boehnke/ExomeComparison/WuSpace_2514360.bed');
    }

    # Default ref seq
    unless (defined $self->reference_fasta) {
        $self->reference_fasta(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa');
    }

    ## Set a default reference transcript annotator version ##
#    my $build_id = $self->build_id;
#    my $build = Genome::Model::Build->get(build_id => $build_id);
my $annotation_reference_transcripts;
if ($self->reference_fasta =~ m/build36/) {
    $annotation_reference_transcripts = "NCBI-human.combined-annotation/54_36p_v3";
}
else {
    $annotation_reference_transcripts = "NCBI-human.combined-annotation/58_37c_v2";
}
#    my $ref_id = $build->model->reference_sequence_build->id;

#    ## Human build 37 ##    
#    if($ref_id == 106942997 || $ref_id == 102671028)
#    {
#        $annotation_reference_transcripts = "NCBI-human.combined-annotation/57_37b";
#    }
#    elsif($ref_id == 104420993 || $ref_id == 103785428)
#    {
#        $annotation_reference_transcripts = "NCBI-mouse.combined-annotation/54_37g_v2";
#    }


    # Set (hardcoded) defaults for tools that have defaults that do not agree with somatic pipeline
#annotation options
    unless (defined $self->annotate_no_headers) {
        $self->annotate_no_headers(1);
    }
    unless (defined $self->transcript_annotation_filter) {
        $self->transcript_annotation_filter("top");
    }
    unless (defined $self->annotation_reference_transcripts) {
        $self->annotation_reference_transcripts($annotation_reference_transcripts);
    }
    unless (defined $self->transcript_variant_annotator_version) {
        $self->transcript_variant_annotator_version(Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version);
    }

#tiering options
    unless (defined $self->only_tier_1) {
        $self->only_tier_1(0);
    }
    unless (defined $self->only_tier_1_indel) {
        $self->only_tier_1_indel(0);
    }
    unless (defined $self->skip_roi) {
        $self->skip_roi(0);
    }

#dbsnp lookup-variants settings
    unless (defined $self->report_mode) {
        $self->report_mode("known-only");
    }
    unless (defined $self->append_rs_id) {
        $self->append_rs_id(1);
    }
    unless (defined $self->dbSNP_version) {
        $self->dbSNP_version(130);
    }



#filter false positives
    unless (defined $self->analysis_type) {
        $self->analysis_type("capture");
    }

#maf file creation
    unless (defined $self->project_name) {
        $self->project_name("Germline Project");
    }
    unless (defined $self->center) {
        $self->center("genome.wustl.edu");
    }
    unless (defined $self->build) {
        $self->build("NCBI-human-build36");
    }
    unless (defined $self->sequence_phase) {
        $self->sequence_phase("4");
    }
    unless (defined $self->sequence_source) {
        $self->sequence_source("Capture");
    }
    unless (defined $self->sequencer) {
        $self->sequencer("Illumina_GAIIx_or_Hiseq");
    }

    return 1;
}

sub default_filenames{
    my $self = shift;
   
    my %default_filenames = (        
        ## samtools Adapted Files ##
        adaptor_output_indel                => 'samtools.output.indel.formatted',
        samtools_snp_output_adaptor         => 'samtools.output.snp.adaptor',

        ## Varscan Output Files ##
        varscan_snp_output                  => 'varScan.output.snp',
        varscan_indel_output                => 'varScan.output.indel',

        ## Varscan Adapted Output Files ##
        varscan_adaptor_snp                 => 'varScan.output.snp.formatted',
        varscan_adaptor_indel               => 'varScan.output.indel.formatted',

        ## GATK Output Files
        GATK_indel_output                   => 'GATK.output.indel',
	GATK_indel_formatted_output         => 'GATK.output.indel.formatted',

        ## GATK Unified Genotyper Output Files
	GATK_ug_indel_output                => 'GATK.ug.output.indel',

        ## GATK Unified Genotyper Adapted Output Files
        GATK_ug_adaptor_indel                  => 'GATK.ug.output.indel.adaptor',

        ## Combined samtools+Varscan Output files ##
        merged_snp_output                   => 'merged.germline.snp',            ## Generated from merge-variants of samtools and varScan
        merged_indel_output                 => 'merged.germline.indel',          ## Generated from merge-variants of samtools and varScan and gatk ##

        ## Limit to ROI, Combined samtools+Varscan Output files ##
        merged_snp_output_ROI               => 'merged.germline.snp.ROI',          ## Generated from merge-variants of samtools and varScan ##
        merged_indel_output_ROI             => 'merged.germline.indel.ROI',          ## Generated from merge-variants of samtools and varScan and gatk ##

        ## Filtered SNP and Indel files ##
        tier_1_snpfilter_file               => 'merged.germline.snp.ROI.strandfilter',
        tier_1_snpfilter_file_filtered      => 'merged.germline.snp.ROI.strandfilter_filtered',
        tier_1_indelfilter_file             => 'merged.germline.indel.ROI.strandfilter',
        tier_1_indelfilter_file_filtered    => 'merged.germline.indel.ROI.strandfilter_filtered',

        ## Annotation output files ##
        annotate_output_snp                 => 'annotation.germline.snp.strandfilter.transcript',
        annotate_output_indel               => 'annotation.germline.indel.strandfilter.transcript',
        ucsc_output_snp                     => 'annotation.germline.snp.strandfilter.ucsc',
        ucsc_output_indel                   => 'annotation.germline.indel.strandfilter.ucsc',
        ucsc_unannotated_output             => 'annotation.germline.snp.strandfilter.unannot-ucsc',
        ucsc_unannotated_indel_output       => 'annotation.germline.indel.strandfilter.unannot-ucsc',

        ## Tiered SNP files (all confidence) ##
        tier_1_snp_file                     => 'merged.germline.snp.ROI.strandfilter.tier1.out',
        tier_2_snp_file                     => 'merged.germline.snp.ROI.strandfilter.tier2.out',
        tier_3_snp_file                     => 'merged.germline.snp.ROI.strandfilter.tier3.out',
        tier_4_snp_file                     => 'merged.germline.snp.ROI.strandfilter.tier4.out',

        ## Tiered indel files (all confidence) ##
        tier_1_indel_file                   => 'merged.germline.indel.ROI.strandfilter.tier1.out',
        tier_2_indel_file                   => 'merged.germline.indel.ROI.strandfilter.tier2.out',
        tier_3_indel_file                   => 'merged.germline.indel.ROI.strandfilter.tier3.out',
        tier_4_indel_file                   => 'merged.germline.indel.ROI.strandfilter.tier4.out',

        ## dbsnp  ##
        tier_1_dbsnp_file                   => 'merged.germline.snp.ROI.tier1.out.dbsnp',

	## maf file ##
        tier1_maf_file                      => 'merged.germline.ROI.tier1.out.maf',

	## vcf file ##
        tier1_vcf_file                      => 'merged.germline.ROI.tier1.out.vcf',

    );

    return %default_filenames;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Germline Pipeline" logDir="/gsc/var/log/genome/germline_capture_pipeline">

<!-- VARSCAN2 GERMLINE -->

<!--  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Varscan Germline" toProperty="skip_if_output_present" /> -->
  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="Varscan Germline" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="reference_fasta" toOperation="Varscan Germline" toProperty="reference" />
  <link fromOperation="input connector" fromProperty="varscan_snp_output" toOperation="Varscan Germline" toProperty="output_snp" />
  <link fromOperation="input connector" fromProperty="varscan_indel_output" toOperation="Varscan Germline" toProperty="output_indel" />

<!-- FORMAT VARSCAN SNPS -->

  <link fromOperation="Varscan Germline" fromProperty="output_snp" toOperation="Format Varscan Snvs" toProperty="variants_file" />
  <link fromOperation="input connector" fromProperty="varscan_adaptor_snp" toOperation="Format Varscan Snvs" toProperty="output_file" />

<!-- FORMAT FILTERED SAMTOOLS SNPS -->

  <link fromOperation="input connector" fromProperty="filtered_indelpe_snps" toOperation="Format Samtools Snvs" toProperty="variants_file" />
  <link fromOperation="input connector" fromProperty="samtools_snp_output_adaptor" toOperation="Format Samtools Snvs" toProperty="output_file" />

<!-- MERGE FILTERED SAMTOOLS SNPS AND VARSCAN SNPS -->

  <link fromOperation="Format Varscan Snvs" fromProperty="output_file" toOperation="Merge SNPs" toProperty="varscan_file" />
  <link fromOperation="Format Samtools Snvs" fromProperty="output_file" toOperation="Merge SNPs" toProperty="glf_file" />
  <link fromOperation="input connector" fromProperty="merged_snp_output" toOperation="Merge SNPs" toProperty="output_file" />

<!-- Limit SNPs ROI --> 

  <link fromOperation="Merge SNPs" fromProperty="output_file" toOperation="Limit SNPs ROI" toProperty="input_file" />
  <link fromOperation="input connector" fromProperty="regions_file" toOperation="Limit SNPs ROI" toProperty="regions_file" />
  <link fromOperation="input connector" fromProperty="merged_snp_output_ROI" toOperation="Limit SNPs ROI" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="skip_roi" toOperation="Limit SNPs ROI" toProperty="skip_roi" />

<!-- FILTER SNP VARIANTS -->
  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Filter Snp" toProperty="skip_if_output_present" />
  <link fromOperation="Limit SNPs ROI" fromProperty="output_file" toOperation="Filter Snp" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="Filter Snp" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="tier_1_snpfilter_file" toOperation="Filter Snp" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="tier_1_snpfilter_file_filtered" toOperation="Filter Snp" toProperty="filtered_file" />
  <link fromOperation="input connector" fromProperty="reference_fasta" toOperation="Filter Snp" toProperty="reference" />

<!-- RUN TRANSCRIPT ANNOTATION FOR SNPS --> 

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Annotate Transcript Variants Snp" toProperty="skip_if_output_present" />
  <link fromOperation="Filter Snp" fromProperty="output_file" toOperation="Annotate Transcript Variants Snp" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="annotate_output_snp" toOperation="Annotate Transcript Variants Snp" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="annotate_no_headers" toOperation="Annotate Transcript Variants Snp" toProperty="no_headers" />
  <link fromOperation="input connector" fromProperty="transcript_annotation_filter" toOperation="Annotate Transcript Variants Snp" toProperty="annotation_filter" />
  <link fromOperation="input connector" fromProperty="annotation_reference_transcripts" toOperation="Annotate Transcript Variants Snp" toProperty="reference_transcripts" />
  <link fromOperation="input connector" fromProperty="transcript_variant_annotator_version" toOperation="Annotate Transcript Variants Snp" toProperty="use_version" />

<!-- RUN UCSC ANNOTATION FOR SNPS --> 

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Annotate UCSC" toProperty="skip_if_output_present" />
  <link fromOperation="Filter Snp" fromProperty="output_file" toOperation="Annotate UCSC" toProperty="input_file" />
  <link fromOperation="input connector" fromProperty="ucsc_output_snp" toOperation="Annotate UCSC" toProperty="output_file" /> 
  <link fromOperation="input connector" fromProperty="ucsc_unannotated_output" toOperation="Annotate UCSC" toProperty="unannotated_file" /> 
  <link fromOperation="input connector" fromProperty="only_tier_1" toOperation="Annotate UCSC" toProperty="skip" /> 

<!-- TIER VARIANTS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Tier Variants Snp" toProperty="skip_if_output_present" />
  <link fromOperation="input connector" fromProperty="tier_1_snp_file" toOperation="Tier Variants Snp" toProperty="tier1_file" />
  <link fromOperation="input connector" fromProperty="tier_2_snp_file" toOperation="Tier Variants Snp" toProperty="tier2_file" />
  <link fromOperation="input connector" fromProperty="tier_3_snp_file" toOperation="Tier Variants Snp" toProperty="tier3_file" />
  <link fromOperation="input connector" fromProperty="tier_4_snp_file" toOperation="Tier Variants Snp" toProperty="tier4_file" />
  <link fromOperation="input connector" fromProperty="only_tier_1" toOperation="Tier Variants Snp" toProperty="only_tier_1" />
  <link fromOperation="Annotate UCSC" fromProperty="output_file" toOperation="Tier Variants Snp" toProperty="ucsc_file" />
  <link fromOperation="Merge SNPs" fromProperty="output_file" toOperation="Tier Variants Snp" toProperty="variant_file" />
  <link fromOperation="Annotate Transcript Variants Snp" fromProperty="output_file" toOperation="Tier Variants Snp" toProperty="transcript_annotation_file" />

<!-- MARK DBSNP VARIANTS -->

  <link fromOperation="Tier Variants Snp" fromProperty="tier1_file" toOperation="dbSNP Snp" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="tier_1_dbsnp_file" toOperation="dbSNP Snp" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="report_mode" toOperation="dbSNP Snp" toProperty="report_mode" />
  <link fromOperation="input connector" fromProperty="append_rs_id" toOperation="dbSNP Snp" toProperty="append_rs_id" />
  <link fromOperation="input connector" fromProperty="dbSNP_version" toOperation="dbSNP Snp" toProperty="dbSNP_version" />

<!-- INDELS START HERE -->

<!-- GATK GERMLINE -->

  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="GATK Germline" toProperty="bam_file" /> 
  <link fromOperation="input connector" fromProperty="GATK_indel_output" toOperation="GATK Germline" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="GATK_indel_formatted_output" toOperation="GATK Germline" toProperty="formatted_file" />
  <link fromOperation="input connector" fromProperty="reference_fasta" toOperation="GATK Germline" toProperty="reference" />    
  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="GATK Germline" toProperty="skip_if_output_present" />    

<!-- GATK Unified Genotyper GERMLINE -->
  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="GATK Unified Genotyper Germline" toProperty="bam_file" /> 
  <link fromOperation="input connector" fromProperty="GATK_ug_indel_output" toOperation="GATK Unified Genotyper Germline" toProperty="vcf_output_file" />
  <link fromOperation="input connector" fromProperty="reference_fasta" toOperation="GATK Unified Genotyper Germline" toProperty="reference_fasta" />    
  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="GATK Unified Genotyper Germline" toProperty="skip_if_output_present" />    

<!-- FORMAT GATK Unified Genotyper INDELS -->

  <link fromOperation="GATK Unified Genotyper Germline" fromProperty="vcf_output_file" toOperation="Format GATK Indels" toProperty="variants_file" /> 
  <link fromOperation="input connector" fromProperty="GATK_ug_adaptor_indel" toOperation="Format GATK Indels" toProperty="output_file" /> 

<!-- FORMAT SAMTOOLS INDELS -->

  <link fromOperation="input connector" fromProperty="indels_all_sequences_filtered" toOperation="Format Samtools Indels" toProperty="variants_file" />
  <link fromOperation="input connector" fromProperty="adaptor_output_indel" toOperation="Format Samtools Indels" toProperty="output_file" />

<!-- FORMAT VARSCAN INDELS -->

  <link fromOperation="Varscan Germline" fromProperty="output_indel" toOperation="Format Varscan Indels" toProperty="variants_file" />
  <link fromOperation="input connector" fromProperty="varscan_adaptor_indel" toOperation="Format Varscan Indels" toProperty="output_file" />

<!-- MERGE ADAPTED INDELS FROM SAMTOOLS AND VARSCAN AND GATK-->

  <link fromOperation="Format Varscan Indels" fromProperty="output_file" toOperation="Merge Indels" toProperty="varscan_file" />
  <link fromOperation="Format Samtools Indels" fromProperty="output_file" toOperation="Merge Indels" toProperty="glf_file" />
  <link fromOperation="GATK Germline" fromProperty="formatted_file" toOperation="Merge Indels" toProperty="gatk_file" />
  <link fromOperation="input connector" fromProperty="merged_indel_output" toOperation="Merge Indels" toProperty="output_file" />
  <link fromOperation="Format GATK Indels" fromProperty="output_file" toOperation="Merge Indels" toProperty="gatk_unified_file" />

<!-- Limit Indels ROI -->

  <link fromOperation="Merge Indels" fromProperty="output_file" toOperation="Limit Indels ROI" toProperty="input_file" />
  <link fromOperation="input connector" fromProperty="regions_file" toOperation="Limit Indels ROI" toProperty="regions_file" />
  <link fromOperation="input connector" fromProperty="merged_indel_output_ROI" toOperation="Limit Indels ROI" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="skip_roi" toOperation="Limit Indels ROI" toProperty="skip_roi" />

<!-- FILTER INDEL VARIANTS -->

  <link fromOperation="Limit Indels ROI" fromProperty="output_file" toOperation="Filter Indel" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="Filter Indel" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="tier_1_indelfilter_file" toOperation="Filter Indel" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="tier_1_indelfilter_file_filtered" toOperation="Filter Indel" toProperty="filtered_file" />
  <link fromOperation="input connector" fromProperty="reference_fasta" toOperation="Filter Indel" toProperty="reference" />

<!-- RUN TRANSCRIPT ANNOTATION FOR INDELS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Annotate Transcript Variants Indel" toProperty="skip_if_output_present" />
  <link fromOperation="Filter Indel" fromProperty="output_file" toOperation="Annotate Transcript Variants Indel" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="annotate_output_indel" toOperation="Annotate Transcript Variants Indel" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="annotate_no_headers" toOperation="Annotate Transcript Variants Indel" toProperty="no_headers" />
  <link fromOperation="input connector" fromProperty="transcript_annotation_filter" toOperation="Annotate Transcript Variants Indel" toProperty="annotation_filter" />

<!-- RUN UCSC ANNOTATION FOR INDELS --> 

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Annotate UCSC Indel" toProperty="skip_if_output_present" />
  <link fromOperation="Filter Indel" fromProperty="output_file" toOperation="Annotate UCSC Indel" toProperty="input_file" />
  <link fromOperation="input connector" fromProperty="ucsc_output_indel" toOperation="Annotate UCSC Indel" toProperty="output_file" /> 
  <link fromOperation="input connector" fromProperty="ucsc_unannotated_indel_output" toOperation="Annotate UCSC Indel" toProperty="unannotated_file" /> 
  <link fromOperation="input connector" fromProperty="only_tier_1_indel" toOperation="Annotate UCSC Indel" toProperty="skip" /> 

<!-- TIER INDELS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Tier Variants Indel" toProperty="skip_if_output_present" />
  <link fromOperation="input connector" fromProperty="tier_1_indel_file" toOperation="Tier Variants Indel" toProperty="tier1_file" />
  <link fromOperation="input connector" fromProperty="tier_2_indel_file" toOperation="Tier Variants Indel" toProperty="tier2_file" />
  <link fromOperation="input connector" fromProperty="tier_3_indel_file" toOperation="Tier Variants Indel" toProperty="tier3_file" />
  <link fromOperation="input connector" fromProperty="tier_4_indel_file" toOperation="Tier Variants Indel" toProperty="tier4_file" />
  <link fromOperation="input connector" fromProperty="only_tier_1_indel" toOperation="Tier Variants Indel" toProperty="only_tier_1" />
  <link fromOperation="Annotate UCSC Indel" fromProperty="output_file" toOperation="Tier Variants Indel" toProperty="ucsc_file" />
  <link fromOperation="Merge Indels" fromProperty="output_file" toOperation="Tier Variants Indel" toProperty="variant_file" />
  <link fromOperation="Annotate Transcript Variants Indel" fromProperty="output_file" toOperation="Tier Variants Indel" toProperty="transcript_annotation_file" />

<!-- SNV AND INDEL TO MAF -->

  <link fromOperation="Tier Variants Snp" fromProperty="tier1_file" toOperation="Tier1 Maf" toProperty="variant_file" />
  <link fromOperation="dbSNP Snp" fromProperty="output_file" toOperation="Tier1 Maf" toProperty="dbsnp_file" />
  <link fromOperation="Filter Snp" fromProperty="output_file" toOperation="Tier1 Maf" toProperty="snv_filtered_file" />
  <link fromOperation="Filter Snp" fromProperty="filtered_file" toOperation="Tier1 Maf" toProperty="snv_failfiltered_file" />
  <link fromOperation="Annotate Transcript Variants Snp" fromProperty="output_file" toOperation="Tier1 Maf" toProperty="snv_annotation_file" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier1_file" toOperation="Tier1 Maf" toProperty="indel_file" />
  <link fromOperation="Filter Indel" fromProperty="output_file" toOperation="Tier1 Maf" toProperty="indel_filtered_file" />
  <link fromOperation="Filter Indel" fromProperty="filtered_file" toOperation="Tier1 Maf" toProperty="indel_failfiltered_file" />
  <link fromOperation="Annotate Transcript Variants Indel" fromProperty="output_file" toOperation="Tier1 Maf" toProperty="indel_annotation_file" />

  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="Tier1 Maf" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="build_id" toOperation="Tier1 Maf" toProperty="build_id" />  
  <link fromOperation="input connector" fromProperty="tier1_maf_file" toOperation="Tier1 Maf" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="project_name" toOperation="Tier1 Maf" toProperty="project_name" />
  <link fromOperation="input connector" fromProperty="center" toOperation="Tier1 Maf" toProperty="center" />
  <link fromOperation="input connector" fromProperty="build" toOperation="Tier1 Maf" toProperty="build" />
  <link fromOperation="input connector" fromProperty="sequence_phase" toOperation="Tier1 Maf" toProperty="sequence_phase" />
  <link fromOperation="input connector" fromProperty="sequence_source" toOperation="Tier1 Maf" toProperty="sequence_source" />
  <link fromOperation="input connector" fromProperty="sequencer" toOperation="Tier1 Maf" toProperty="sequencer" />

<!-- SNV AND INDEL TO VCF -->

  <link fromOperation="Tier Variants Snp" fromProperty="tier1_file" toOperation="Tier1 Vcf" toProperty="variant_file" />
  <link fromOperation="dbSNP Snp" fromProperty="output_file" toOperation="Tier1 Vcf" toProperty="dbsnp_file" />
  <link fromOperation="Filter Snp" fromProperty="output_file" toOperation="Tier1 Vcf" toProperty="snv_filtered_file" />
  <link fromOperation="Filter Snp" fromProperty="filtered_file" toOperation="Tier1 Vcf" toProperty="snv_failfiltered_file" />
  <link fromOperation="Annotate Transcript Variants Snp" fromProperty="output_file" toOperation="Tier1 Vcf" toProperty="snv_annotation_file" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier1_file" toOperation="Tier1 Vcf" toProperty="indel_file" />
  <link fromOperation="Filter Indel" fromProperty="output_file" toOperation="Tier1 Vcf" toProperty="indel_filtered_file" />
  <link fromOperation="Filter Indel" fromProperty="filtered_file" toOperation="Tier1 Vcf" toProperty="indel_failfiltered_file" />
  <link fromOperation="Annotate Transcript Variants Indel" fromProperty="output_file" toOperation="Tier1 Vcf" toProperty="indel_annotation_file" />

  <link fromOperation="input connector" fromProperty="filtered_indelpe_snps" toOperation="Tier1 Vcf" toProperty="samtools_file" />
  <link fromOperation="Varscan Germline" fromProperty="output_snp" toOperation="Tier1 Vcf" toProperty="varscan_file" />
  <link fromOperation="input connector" fromProperty="germline_bam_file" toOperation="Tier1 Vcf" toProperty="bam_file" />
  <link fromOperation="input connector" fromProperty="build_id" toOperation="Tier1 Vcf" toProperty="build_id" />  
  <link fromOperation="input connector" fromProperty="tier1_vcf_file" toOperation="Tier1 Vcf" toProperty="output_file" />
  <link fromOperation="input connector" fromProperty="project_name" toOperation="Tier1 Vcf" toProperty="project_name" />
  <link fromOperation="input connector" fromProperty="center" toOperation="Tier1 Vcf" toProperty="center" />
  <link fromOperation="input connector" fromProperty="build" toOperation="Tier1 Vcf" toProperty="build" />
  <link fromOperation="input connector" fromProperty="sequence_phase" toOperation="Tier1 Vcf" toProperty="sequence_phase" />
  <link fromOperation="input connector" fromProperty="sequence_source" toOperation="Tier1 Vcf" toProperty="sequence_source" />
  <link fromOperation="input connector" fromProperty="sequencer" toOperation="Tier1 Vcf" toProperty="sequencer" />

<!-- OUTPUT CONNECTORS -->

  <link fromOperation="Tier1 Maf" fromProperty="output_file" toOperation="output connector" toProperty="maf_file_out" />
  <link fromOperation="Tier1 Vcf" fromProperty="output_file" toOperation="output connector" toProperty="vcf_file_out" />

  <link fromOperation="Tier Variants Snp" fromProperty="tier2_file" toOperation="output connector" toProperty="tier_2_snp" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier3_file" toOperation="output connector" toProperty="tier_3_snp" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier4_file" toOperation="output connector" toProperty="tier_4_snp" />

  <link fromOperation="Tier Variants Indel" fromProperty="tier2_file" toOperation="output connector" toProperty="tier_2_indel_output" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier3_file" toOperation="output connector" toProperty="tier_3_indel_output" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier4_file" toOperation="output connector" toProperty="tier_4_indel_output" />

  <operation name="Varscan Germline">
    <operationtype commandClass="Genome::Model::Tools::Varscan::Germline" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Format Varscan Snvs">
    <operationtype commandClass="Genome::Model::Tools::Capture::FormatSnvs" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Format Samtools Snvs">
    <operationtype commandClass="Genome::Model::Tools::Capture::FormatSnvs" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Merge SNPs">
    <operationtype commandClass="Genome::Model::Tools::Capture::MergeVariantCalls" typeClass="Workflow::OperationType::Command" />
  </operation>  

  <operation name="GATK Germline">
    <operationtype commandClass="Genome::Model::Tools::Gatk::GermlineIndel" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="GATK Unified Genotyper Germline">
    <operationtype commandClass="Genome::Model::Tools::Gatk::GermlineIndelUnifiedGenotyper" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Format Varscan Indels">
    <operationtype commandClass="Genome::Model::Tools::Capture::FormatIndels" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Format Samtools Indels">
    <operationtype commandClass="Genome::Model::Tools::Capture::FormatIndels" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Format GATK Indels">
    <operationtype commandClass="Genome::Model::Tools::Gatk::FormatVcf" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Merge Indels">
    <operationtype commandClass="Genome::Model::Tools::Capture::MergeAdaptedIndels" typeClass="Workflow::OperationType::Command" />
  </operation>  

  <operation name="Limit SNPs ROI">
    <operationtype commandClass="Genome::Model::Tools::Capture::LimitToRoi" typeClass="Workflow::OperationType::Command" />
  </operation> 
  <operation name="Limit Indels ROI">
    <operationtype commandClass="Genome::Model::Tools::Capture::LimitToRoi" typeClass="Workflow::OperationType::Command" />
  </operation> 

  <operation name="Annotate UCSC">
      <operationtype commandClass="Genome::Model::Tools::Somatic::UcscAnnotator" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Annotate UCSC Indel">
      <operationtype commandClass="Genome::Model::Tools::Somatic::UcscAnnotator" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Annotate Transcript Variants Snp">
    <operationtype commandClass="Genome::Model::Tools::Annotate::TranscriptVariants" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Annotate Transcript Variants Indel">
    <operationtype commandClass="Genome::Model::Tools::Annotate::TranscriptVariants" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Tier Variants Snp">
    <operationtype commandClass="Genome::Model::Tools::Somatic::TierVariants" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Tier Variants Indel">
    <operationtype commandClass="Genome::Model::Tools::Somatic::TierVariants" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="dbSNP Snp">
    <operationtype commandClass="Genome::Model::Tools::Annotate::LookupVariants" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Filter Snp">
    <operationtype commandClass="Genome::Model::Tools::Somatic::FilterFalsePositives" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Filter Indel">
    <operationtype commandClass="Genome::Model::Tools::Somatic::FilterFalseIndels" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Tier1 Maf">
    <operationtype commandClass="Genome::Model::Tools::Germline::MafMaker" typeClass="Workflow::OperationType::Command" />
  </operation>
  <operation name="Tier1 Vcf">
    <operationtype commandClass="Genome::Model::Tools::Germline::VcfMaker" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>build_id</inputproperty>
    <inputproperty>filtered_indelpe_snps</inputproperty>
    <inputproperty>indels_all_sequences_filtered</inputproperty>
    <inputproperty>regions_file</inputproperty>

    <inputproperty isOptional="Y">germline_bam_file</inputproperty>
    <inputproperty isOptional="Y">reference_fasta</inputproperty>

    <inputproperty isOptional="Y">samtools_snp_output_adaptor</inputproperty>
    <inputproperty isOptional="Y">adaptor_output_indel</inputproperty>

    <inputproperty isOptional="Y">skip_if_output_present</inputproperty>

    <inputproperty isOptional="Y">ucsc_output_snp</inputproperty>
    <inputproperty isOptional="Y">ucsc_output_indel</inputproperty>
    
    <inputproperty isOptional="Y">only_tier_1</inputproperty>
    <inputproperty isOptional="Y">only_tier_1_indel</inputproperty>

    <inputproperty isOptional="Y">data_directory</inputproperty>

    <inputproperty isOptional="Y">annotate_output_germline_snp</inputproperty>
    <inputproperty isOptional="Y">annotate_output_germline_indel</inputproperty>

    <inputproperty isOptional="Y">annotate_output_indel</inputproperty>
    <inputproperty isOptional="Y">annotate_output_snp</inputproperty>
    <inputproperty isOptional="Y">annotate_no_headers</inputproperty>
    <inputproperty isOptional="Y">transcript_annotation_filter</inputproperty>
    <inputproperty isOptional="Y">annotation_reference_transcripts</inputproperty>
    <inputproperty isOptional="Y">transcript_variant_annotator_version</inputproperty>

    <inputproperty isOptional="Y">ucsc_unannotated_output</inputproperty>
    <inputproperty isOptional="Y">ucsc_unannotated_indel_output</inputproperty>

    <inputproperty isOptional="Y">varscan_snp_output</inputproperty>
    <inputproperty isOptional="Y">varscan_snp_germline</inputproperty>
    <inputproperty isOptional="Y">varscan_snp_loh</inputproperty>
    <inputproperty isOptional="Y">varscan_snp_germline</inputproperty>
    <inputproperty isOptional="Y">varscan_indel_output</inputproperty>
    <inputproperty isOptional="Y">varscan_indel_germline</inputproperty>
    <inputproperty isOptional="Y">varscan_indel_loh</inputproperty>
    <inputproperty isOptional="Y">varscan_indel_germline</inputproperty>
    <inputproperty isOptional="Y">varscan_adaptor_snp</inputproperty>
    <inputproperty isOptional="Y">varscan_adaptor_indel</inputproperty>

    <inputproperty isOptional="Y">GATK_indel_output</inputproperty>
    <inputproperty isOptional="Y">GATK_indel_formatted_output</inputproperty>
    <inputproperty isOptional="Y">GATK_ug_indel_output</inputproperty>
    <inputproperty isOptional="Y">GATK_ug_adaptor_indel</inputproperty>

    <inputproperty isOptional="Y">merged_snp_output</inputproperty>
    <inputproperty isOptional="Y">merged_indel_output</inputproperty>
    <inputproperty isOptional="Y">merged_snp_output_ROI</inputproperty>
    <inputproperty isOptional="Y">merged_indel_output_ROI</inputproperty>
    <inputproperty isOptional="Y">skip_roi</inputproperty>

    <inputproperty isOptional="Y">tier_1_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_2_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_3_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_4_snp_file</inputproperty>
    
    <inputproperty isOptional="Y">tier_1_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_2_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_3_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_4_indel_file</inputproperty>

    <inputproperty isOptional="Y">tier_1_dbsnp_file</inputproperty>
    <inputproperty isOptional="Y">report_mode</inputproperty>
    <inputproperty isOptional="Y">append_rs_id</inputproperty>
    <inputproperty isOptional="Y">dbSNP_version</inputproperty>

    <inputproperty isOptional="Y">tier_1_snpfilter_file</inputproperty>
    <inputproperty isOptional="Y">tier_1_snpfilter_file_filtered</inputproperty>
    <inputproperty isOptional="Y">analysis_type</inputproperty>

    <inputproperty isOptional="Y">tier_1_indelfilter_file</inputproperty>
    <inputproperty isOptional="Y">tier_1_indelfilter_file_filtered</inputproperty>

    <inputproperty isOptional="Y">tier1_maf_file</inputproperty>
    <inputproperty isOptional="Y">tier1_vcf_file</inputproperty>
    <inputproperty isOptional="Y">project_name</inputproperty>
    <inputproperty isOptional="Y">center</inputproperty>
    <inputproperty isOptional="Y">build</inputproperty>
    <inputproperty isOptional="Y">sequence_phase</inputproperty>
    <inputproperty isOptional="Y">sequence_source</inputproperty>
    <inputproperty isOptional="Y">sequencer</inputproperty>

    <outputproperty>maf_file_out</outputproperty>
    <outputproperty>vcf_file_out</outputproperty>

    <outputproperty>tier_2_snp</outputproperty>
    <outputproperty>tier_3_snp</outputproperty>
    <outputproperty>tier_4_snp</outputproperty>

    <outputproperty>tier_2_indel_output</outputproperty>
    <outputproperty>tier_3_indel_output</outputproperty>
    <outputproperty>tier_4_indel_output</outputproperty>

  </operationtype>

</workflow>


