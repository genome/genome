package Genome::Model::Tools::Germline::TierVariants;

########################################################################################################################
# TierVariants.pm - A module for tiering variants
#
#   TO-DO LIST
#   -Run Sample QC checks (SNP and CNV)
#   
########################################################################################################################

use strict;
use warnings;
use Workflow;

class Genome::Model::Tools::Germline::TierVariants {
    is => ['Workflow::Operation::Command'],
    workflow => sub { Workflow::Operation->create_from_xml(\*DATA); }
};

sub help_brief {
    "Runs the **variant tiering, circos, and report** pipeline workflow."
}

sub help_synopsis{
    my $self = shift;
    return <<"EOS"

example:
gmt germline tier-variants --annotate-output-germline-indel=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.indels.annotation --build-id=0 --annotate-output-germline-snp=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.annotation --data-directory=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092_tiering --ucsc-output-snp=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.annotation.ucsc --filtered-indelpe-snps=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.annotation --adapted-indel-file=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.filtered.indels --ucsc-output-indel=/gscmnt/sata424/info/medseq/Freimer-Boehnke/FB_2_100_dedup/H_HY-02092.indels.annotation.ucsc

EOS
}

sub help_detail {
    my $self = shift;
    return <<"EOS"
This tool runs a pipeline to take in annotated SNPs and indels. It results in tiered SNPs and indels.
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

    # Set (hardcoded) defaults for tools that have defaults that do not agree with somatic pipeline
    unless (defined $self->skip_if_output_present) {
        $self->skip_if_output_present(1);
    }
    unless (defined $self->only_tier_1) {
        $self->only_tier_1(0);
    }
    unless (defined $self->only_tier_1_indel) {
        $self->only_tier_1_indel(0);
    }

    return 1;
}

sub default_filenames{
    my $self = shift;
   
    my %default_filenames = (        
        ## Tiered SNP and indel files (all confidence) ##
        tier_1_snp_file                     => 'merged.germline.snp.tier1.out',
        tier_2_snp_file                     => 'merged.germline.snp.tier2.out',
        tier_3_snp_file                     => 'merged.germline.snp.tier3.out',
        tier_4_snp_file                     => 'merged.germline.snp.tier4.out',
        tier_1_indel_file                   => 'merged.germline.indel.tier1.out',
        tier_2_indel_file                   => 'merged.germline.indel.tier2.out',
        tier_3_indel_file                   => 'merged.germline.indel.tier3.out',
        tier_4_indel_file                   => 'merged.germline.indel.tier4.out',

        ## Other pipeline output files ##
        circos_graph                        => 'circos_graph.out',
        variant_report_output               => 'cancer_report.html', 
        file_summary_report_output          => 'file_summary_report.html',
    );

    return %default_filenames;
}

1;
__DATA__
<?xml version='1.0' standalone='yes'?>

<workflow name="Tiering Pipeline" logDir="/gsc/var/log/genome/tiering_pipeline">

<!-- TIER VARIANTS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Tier Variants Snp" toProperty="skip_if_output_present" />
  <link fromOperation="input connector" fromProperty="tier_1_snp_file" toOperation="Tier Variants Snp" toProperty="tier1_file" />
  <link fromOperation="input connector" fromProperty="tier_2_snp_file" toOperation="Tier Variants Snp" toProperty="tier2_file" />
  <link fromOperation="input connector" fromProperty="tier_3_snp_file" toOperation="Tier Variants Snp" toProperty="tier3_file" />
  <link fromOperation="input connector" fromProperty="tier_4_snp_file" toOperation="Tier Variants Snp" toProperty="tier4_file" />
  <link fromOperation="input connector" fromProperty="only_tier_1" toOperation="Tier Variants Snp" toProperty="only_tier_1" />
  <link fromOperation="input connector" fromProperty="ucsc_output_snp" toOperation="Tier Variants Snp" toProperty="ucsc_file" />
  <link fromOperation="input connector" fromProperty="filtered_indelpe_snps" toOperation="Tier Variants Snp" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="annotate_output_germline_snp" toOperation="Tier Variants Snp" toProperty="transcript_annotation_file" />

<!-- TIER INDELS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Tier Variants Indel" toProperty="skip_if_output_present" />
  <link fromOperation="input connector" fromProperty="tier_1_indel_file" toOperation="Tier Variants Indel" toProperty="tier1_file" />
  <link fromOperation="input connector" fromProperty="tier_2_indel_file" toOperation="Tier Variants Indel" toProperty="tier2_file" />
  <link fromOperation="input connector" fromProperty="tier_3_indel_file" toOperation="Tier Variants Indel" toProperty="tier3_file" />
  <link fromOperation="input connector" fromProperty="tier_4_indel_file" toOperation="Tier Variants Indel" toProperty="tier4_file" />
  <link fromOperation="input connector" fromProperty="only_tier_1_indel" toOperation="Tier Variants Indel" toProperty="only_tier_1" />
  <link fromOperation="input connector" fromProperty="ucsc_output_indel" toOperation="Tier Variants Indel" toProperty="ucsc_file" />
  <link fromOperation="input connector" fromProperty="adapted_indel_file" toOperation="Tier Variants Indel" toProperty="variant_file" />
  <link fromOperation="input connector" fromProperty="annotate_output_germline_indel" toOperation="Tier Variants Indel" toProperty="transcript_annotation_file" />

<!-- PLOT CIRCOS -->

  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Plot Circos" toProperty="skip_if_output_present" />
  <link fromOperation="input connector" fromProperty="circos_graph" toOperation="Plot Circos" toProperty="output_file" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier1_file" toOperation="Plot Circos" toProperty="tier1_hclabel_file" />

<!-- WAIT FOR CIRCOS -->

  <link fromOperation="input connector" fromProperty="build_id" toOperation="Wait for Circos" toProperty="build_id" />
  <link fromOperation="Plot Circos" fromProperty="result" toOperation="Wait for Circos" toProperty="plot circos result" />

<!-- GENERATE REPORT -->
 
  <link fromOperation="Wait for Circos" fromProperty="build_id" toOperation="Generate Reports" toProperty="build_id" />
  <link fromOperation="input connector" fromProperty="variant_report_output" toOperation="Generate Reports" toProperty="variant_report_output" />
  <link fromOperation="input connector" fromProperty="file_summary_report_output" toOperation="Generate Reports" toProperty="file_summary_report_output" />
  <link fromOperation="input connector" fromProperty="skip_if_output_present" toOperation="Generate Reports" toProperty="skip_if_output_present" />

<!-- OUTPUT CONNECTORS -->

  <link fromOperation="Plot Circos" fromProperty="output_file" toOperation="output connector" toProperty="circos_big_graph" />
  <link fromOperation="Generate Reports" fromProperty="variant_report_output" toOperation="output connector" toProperty="final_variant_report_output" />

  <link fromOperation="Tier Variants Snp" fromProperty="tier1_file" toOperation="output connector" toProperty="tier_1_snp" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier2_file" toOperation="output connector" toProperty="tier_2_snp" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier3_file" toOperation="output connector" toProperty="tier_3_snp" />
  <link fromOperation="Tier Variants Snp" fromProperty="tier4_file" toOperation="output connector" toProperty="tier_4_snp" />

  <link fromOperation="Tier Variants Indel" fromProperty="tier1_file" toOperation="output connector" toProperty="tier_1_indel_output" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier2_file" toOperation="output connector" toProperty="tier_2_indel_output" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier3_file" toOperation="output connector" toProperty="tier_3_indel_output" />
  <link fromOperation="Tier Variants Indel" fromProperty="tier4_file" toOperation="output connector" toProperty="tier_4_indel_output" />

  <operation name="Tier Variants Snp">
    <operationtype commandClass="Genome::Model::Tools::Somatic::TierVariants" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Tier Variants Indel">
    <operationtype commandClass="Genome::Model::Tools::Somatic::TierVariants" typeClass="Workflow::OperationType::Command" />
  </operation>


  <operation name="Plot Circos">
    <operationtype commandClass="Genome::Model::Tools::Somatic::PlotCircos" typeClass="Workflow::OperationType::Command" />
  </operation>

  <operation name="Wait for Circos">
      <operationtype typeClass="Workflow::OperationType::Block">
        <property>build_id</property>
        <property>plot circos result</property>
    </operationtype>
  </operation>

  <operation name="Generate Reports">
    <operationtype commandClass="Genome::Model::Tools::Somatic::RunReports" typeClass="Workflow::OperationType::Command" />
  </operation>  

  <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>build_id</inputproperty>
    <inputproperty>filtered_indelpe_snps</inputproperty>
    <inputproperty>adapted_indel_file</inputproperty>
    <inputproperty>ucsc_output_snp</inputproperty>
    <inputproperty>ucsc_output_indel</inputproperty>
    
    <inputproperty isOptional="Y">skip_if_output_present</inputproperty>

    <inputproperty isOptional="Y">only_tier_1</inputproperty>
    <inputproperty isOptional="Y">only_tier_1_indel</inputproperty>

    <inputproperty isOptional="Y">data_directory</inputproperty>

    <inputproperty>annotate_output_germline_snp</inputproperty>
    <inputproperty>annotate_output_germline_indel</inputproperty>



    <inputproperty isOptional="Y">tier_1_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_2_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_3_snp_file</inputproperty>
    <inputproperty isOptional="Y">tier_4_snp_file</inputproperty>
    
    <inputproperty isOptional="Y">tier_1_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_2_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_3_indel_file</inputproperty>
    <inputproperty isOptional="Y">tier_4_indel_file</inputproperty>
   
    <inputproperty isOptional="Y">circos_graph</inputproperty>

    <inputproperty isOptional="Y">variant_report_output</inputproperty>
    <inputproperty isOptional="Y">file_summary_report_output</inputproperty>

    <outputproperty>tier_1_snp</outputproperty>
    <outputproperty>tier_2_snp</outputproperty>
    <outputproperty>tier_3_snp</outputproperty>
    <outputproperty>tier_4_snp</outputproperty>

    <outputproperty>tier_1_indel_output</outputproperty>
    <outputproperty>tier_2_indel_output</outputproperty>
    <outputproperty>tier_3_indel_output</outputproperty>
    <outputproperty>tier_4_indel_output</outputproperty>
    <outputproperty>circos_big_graph</outputproperty>
    <outputproperty>final_variant_report_output</outputproperty>
  </operationtype>

</workflow>


