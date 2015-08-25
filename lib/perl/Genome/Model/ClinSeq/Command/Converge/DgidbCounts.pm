package Genome::Model::ClinSeq::Command::Converge::DgidbCounts;
#Written by Malachi Griffith & Obi Griffith

use strict;
use warnings;
use Genome;
use File::Spec;

class Genome::Model::ClinSeq::Command::Converge::DgidbCounts {
  is => ['Genome::Model::ClinSeq::Command::Converge::Base',
    'Genome::Model::ClinSeq::Util'],
  has_input => [
    outdir => {
      is => 'FilesystemPath',
      doc => 'Directory where output files will be written',
    },
  ],
  doc => 'summarize the number of drug
  interactions found for the following event types',
};

sub help_synopsis {
  return <<INFO;
  Example usage:
genome model clin-seq converge dgidb-counts --builds='id in ["4b7539bb10cc4b9c97577cf11f4c79a2","cdca0edf526c4fe193d3054627a5871b"]' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-counts --builds='model.model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-counts --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge dgidb-counts --builds='model.id in ["279f50e35d2b479ea3c32486eafd4ad4","7143119a93984056ae3f32c88c9ac2a1"],is_last_complete=1' --outdir=/tmp/snv_indel_report/
INFO
}

sub help_detail {
  return <<EOS
  For a group of ClinSeq models, summarize the number of drug interactions found for the following event types:

  1. WGS SNVs, WGS Indels, WGS amplifications
  2. Exome SNVs, Exome Indels
  3. RNA-seq outlier genes

  Longer term, this concept should use the new DGIDB command line modules for getting interactions directly from the database
EOS
}

sub execute {
  my $self = shift;
  my @clinseq_builds = $self->builds;

  #Get files:
  # - drug-gene interaction files for each event type
  # - annotated files for each event type containing potentially druggable results
  my $files = $self->getFiles('-builds'=>\@clinseq_builds);

  #Known druggable genes
  my @event_types = qw(snv indel cnv_gain rna_cufflinks_absolute rna_tophat_absolute);
  my $k_result = $self->parseKnownDruggableFiles('-files'=>$files, '-event_types'=>\@event_types);
  my $k_g = $k_result->{'genes'}; #Known druggable genes
  my $k_i = $k_result->{'interactions'}; #Known druggable gene interactions
  my $data_type_sum = $k_result->{'data_type_sum'};
  my (%totals, %data_type_sum_collapsed);
  $self->write_summary(\%totals, $data_type_sum, \%data_type_sum_collapsed, $k_g);
  $self->write_totals_genes(\%totals);
  return 1;
}

sub write_totals_genes {
  my $self = shift;
  my $totals = shift;
  my $outdir = $self->outdir;
  my $totals_f = File::Spec->catfile(
    $outdir,
    "ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_totals.txt");
  my $genes_f = File::Spec->catfile(
    $outdir,
    "ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary_genes.txt");
  my @totals_header = ("patient", "data_type", "druggable_gene_count");
  my $totals_writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $totals_f,
    separator => "\t",
    headers => \@totals_header,
    print_headers => 1,
  );
  my @genes_header = ("patient", "data_type", "druggable_gene");
  my $genes_writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $genes_f,
    separator => "\t",
    headers => \@genes_header,
    print_headers => 1,
  );
  my $totals_data;
  my $genes_data;
  foreach my $patient (sort keys %$totals){
    foreach my $data_type (sort keys %{$totals->{$patient}}){
      my $gene_count = keys %{$totals->{$patient}{$data_type}};
      $totals_data->{patient} = $patient;
      $totals_data->{data_type} = $data_type;
      $totals_data->{druggable_gene_count} = $gene_count;
      $totals_writer->write_one($totals_data);
      foreach my $gene (sort keys %{$totals->{$patient}{$data_type}}){
        $genes_data->{patient} = $patient;
        $genes_data->{data_type} = $data_type;
        $genes_data->{druggable_gene} = $gene;
        $genes_writer->write_one($genes_data);
      }
    }
  }
}

#Print summary of druggable genes broken down by event and data type
sub write_summary {
  my $self = shift;
  my $totals = shift;
  my $data_type_sum = shift;
  my $data_type_sum_collapsed = shift;
  my $k_g = shift;
  my $outdir = $self->outdir;
  my $outfile = "$outdir/ClinSeq_WGS-EXOME-RNASEQ_drugabble_targets_summary.txt";
  my @header = ("patient", "event_type", "data_type", "gene", "drug");
  my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $outfile,
    separator => "\t",
    headers => \@header,
    print_headers => 1,
  );
  my $out_data;
  foreach my $patient (sort keys %{$data_type_sum}){
    foreach my $event_type (sort keys %{$data_type_sum->{$patient}}){
      foreach my $data_type (sort keys %{$data_type_sum->{$patient}->{$event_type}}){
        foreach my $gene (sort keys %{$data_type_sum->{$patient}->{$event_type}->{$data_type}}){
          #Get drug list for gene
          my %drugs = %{$k_g->{$gene}->{drug_list}};
          foreach my $drug (sort keys %drugs){
            $totals->{$patient}{$data_type}{$gene}++;
            $out_data->{patient} = $patient;
            $out_data->{event_type} = $event_type;
            $out_data->{data_type} = $data_type;
            $out_data->{gene} = $gene;
            $out_data->{drug} = $drug;
            $writer->write_one($out_data);
            my $collapsed_string="$patient\t$event_type\t$data_type\t$gene\t$drug";
            $data_type_sum_collapsed->{$collapsed_string} = 1;
          }
        }
      }
    }
  }
  return $data_type_sum_collapsed;
}

sub get_name {
  my $self = shift;
  my $b = shift;
  my $subject_common_name = "null";
  if($b->model->wgs_model) {
    $subject_common_name = $b->model->wgs_model->last_succeeded_build->subject->common_name;
  } elsif($b->model->exome_model) {
    my $subject_common_name = $b->model->exome_model->last_succeeded_build->subject->common_name;
  }
  return $subject_common_name;
}

1;

