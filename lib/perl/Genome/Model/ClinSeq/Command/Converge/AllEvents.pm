package Genome::Model::ClinSeq::Command::Converge::AllEvents;
use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::AllEvents {
    is => 'Genome::Model::ClinSeq::Command::Converge::Base',
    has_input => [
        outdir => { 
               is => 'FilesystemPath',
               doc => 'Directory where output files will be written',
        },
    ],
    has_optional_input => [
        snv_label => {
                is => 'Text',
                doc => 'Consider SNVs and label affected samples with this string.  e.g. S, SNV, etc.',
        },
        indel_label => {
                is => 'Text',
                doc => 'Consider Indels and label affected samples with this string.  e.g. I, Indel, etc.',
        },
        cnv_gain_label => {
                is => 'Text',
                doc => 'Consider DNA CNV gains (amplifications) and label affected samples with this string.  e.g. A, CnGain, etc.',
        },
        cnv_loss_label => {
                is => 'Text',
                doc => 'Consider DNA CNV losses (deletions) and label affected samples with this string.  e.g. D, CnLoss, etc.',
        },
        rna_de_up => {
                is => 'Text',
                doc => 'Consider RNA DE gains and label affected samples with this string.  e.g. G, RnaUp, etc.',
        },
        rna_de_down => {
                is => 'Text',
                doc => 'Consider RNA DE losses and label affected samples with this string.  e.g. L, RnaUp, etc.',
        },
        target_gene_list => {
                is => 'FilesystemPath',
                doc => 'Limit consideration to genes in this list (Ensembl Gene IDs)',
        },
        ignore_gene_list => {
                is => 'FilesystemPath',
                doc => 'Regardless of anything else do not allow any genes in this list (Ensembl Gene IDs)',
        },
    ],   
    doc => 'converge SNV, InDel, CNV, Exp, and DE results from mutiple clinseq builds into a single table',
};

sub help_synopsis {
  return <<EOS

genome model clin-seq converge all-events --builds='id in [133577030,133611960]' --outdir=/tmp/converge_all_events/

genome model clin-seq converge all-events --builds='model.model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_all_events/

genome model clin-seq converge all-events --builds='model_groups.id=50714,is_last_complete=1' --outdir=/tmp/converge_all_events/

EOS
}

sub execute {
  my $self = shift;
  my @builds = $self->builds;
  my $outdir = $self->outdir;
  $outdir .= "/" unless ($outdir =~ /\/$/);

  #Get human readable names hash, keyed on build id
  my $labels = $self->resolve_clinseq_subject_labels;

  #Produce a table that summarizes the molecular events of each gene, broken down by event type for each clinseq subject
  #Include summary columns that provide samples counts affected by: SNVs, Indels, ..., Any Event Type
  #One gene/transcript per row
  #Each sample column will contain a list of observed event types, each indicated by a single letter, word, etc.

  #For example:
  #S = SNV (non-silent substitution)
  #I = Indel (small indel)
  #A = Amplified (CNV gain)
  #D = Deleted (CNV loss)
  #O = Outlier expression 
  #H = DE up
  #L = DE down

  #The user will specify these as options. Indicating it will trigger that kind of event consideration

  #The final table can optionally be limited to only those genes in a pre-approved list

  #The final table can optionally be filtered to exclude a genes in a user supplied list


  return 1;
};

1;

