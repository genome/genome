package Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport;

use strict;
use warnings;
use Genome;
use Data::Dumper;
use File::Basename;

class Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport {
    is => 'Command::V2',
    has_input => [
        wgs_build           => { is => 'Genome::Model::Build', is_optional => 1, },
        exome_build         => { is => 'Genome::Model::Build', is_optional => 1, },
        tumor_rnaseq_build  => { is => 'Genome::Model::Build', is_optional => 1, },
        normal_rnaseq_build => { is => 'Genome::Model::Build', is_optional => 1, },
       
        cancer_annotation_db => { 
            is => 'Genome::Db::Tgi::CancerAnnotation', 
            example_values => [$Genome::Model::ClinSeq::DEFAULT_CANCER_ANNOTATION_DB_ID],
        },

        wgs_positions_file          => { is => 'FilesystemPath', is_optional => 1 },
        exome_positions_file        => { is => 'FilesystemPath', is_optional => 1 },
        wgs_exome_positions_file    => { is => 'FilesystemPath', is_optional => 1 },
        
        tumor_fpkm_file             => { is => 'FilesystemPath', is_optional => 1 },
    ],
    has_param => [
        verbose => { is => 'Boolean', default_value => 0 },
    ],
    doc => 'Get BAM red counts for SNV positions from WGS, Exome and RNAseq BAMS',     
};

sub positions_files {
    my $self = shift;
    grep { $_ } map { $self->$_ } map { $_ . '_positions_file' } qw/wgs exome wgs_exome/;
}

sub execute {
  my $self = shift;
  $self->status_message("starting summarize tier1 snvs with " . Data::Dumper::Dumper($self));
  my $wgs_build = $self->wgs_build;
  my $exome_build = $self->exome_build;
  my $tumor_rnaseq_build = $self->tumor_rnaseq_build;
  my $normal_rnaseq_build = $self->normal_rnaseq_build;
  my @positions_files = $self->positions_files;
  my $tumor_fpkm_file = $self->tumor_fpkm_file;
  my $verbose = $self->verbose;
  my $cancer_annotation_db = $self->cancer_annotation_db;

  my $read_counts_summary_script = __FILE__ . '.R'; #"$script_dir"."snv/WGS_vs_Exome_vs_RNAseq_VAF_and_FPKM.R";

  $self->status_message("Positions files are " . Data::Dumper::Dumper(\@positions_files));

  foreach my $positions_file (@positions_files){
    my $fb = &getFilePathBase('-path'=>$positions_file);
    my $output_file = $fb->{$positions_file}->{base} . ".readcounts" . $fb->{$positions_file}->{extension};
    my $output_stats_dir = $fb->{$positions_file}->{base_dir} . "summary/";

    my @params = ('positions_file' => $positions_file);
    push (@params, ('wgs_som_var_build' => $wgs_build)) if $wgs_build;
    push (@params, ('exome_som_var_build' => $exome_build)) if $exome_build;
    push (@params, ('rna_seq_tumor_build' => $tumor_rnaseq_build)) if $tumor_rnaseq_build;
    push (@params, ('rna_seq_normal_build' => $normal_rnaseq_build)) if $normal_rnaseq_build;
    push (@params, ('output_file' => $output_file));
    push (@params, ('cancer_annotation_db' => $cancer_annotation_db));
    push (@params, ('verbose' => $verbose));

    $self->status_message("Params for GetBamReadCounts are " . Data::Dumper::Dumper({ @params }));
    my $bam_rc_cmd = Genome::Model::ClinSeq::Command::GetBamReadCounts->create(@params);

    #Summarize the positions file using an R script.  BUT if no variants are present, skip this positions file.
    my $positions_count = 0;
    open (POS, "$positions_file") || die "\n\nCould not open positions file: $positions_file\n\n";
    my $header = 1;
    while(<POS>){
      if ($header){
        $header = 0;
        next();
      }
      $positions_count++;
    }
    close(POS);
    unless($positions_count > 0){
      if ($verbose){
        self->status_message("\n\nNo SNV positions found, skipping summary");
      }
      next();
    }

    #First get the read counts for the current file of SNVs (from WGS, Exome, or WGS+Exome)
    my $r = $bam_rc_cmd->execute();
    unless ($r) {
        $self->error_message("Error from GetBamReadCounts: " . $bam_rc_cmd->error_message);
    }

    #Set up the read count summary script command (an R script)
    my $rc_summary_cmd;
    my $rc_summary_stdout = "$output_stats_dir"."rc_summary.stdout";
    my $rc_summary_stderr = "$output_stats_dir"."rc_summary.stderr";
    if ($tumor_rnaseq_build){
      #my $tumor_fpkm_file = $out_paths->{'tumor_rnaseq_cufflinks_absolute'}->{'isoforms.merged.fpkm.expsort.tsv'}->{path};
      $rc_summary_cmd = "$read_counts_summary_script $output_stats_dir $output_file $tumor_fpkm_file";
    }else{
      $rc_summary_cmd = "$read_counts_summary_script $output_stats_dir $output_file";
    }
    $rc_summary_cmd .= " 1>$rc_summary_stdout 2>$rc_summary_stderr";

    #Summarize the BAM readcounts results for candidate variants - produce descriptive statistics, figures etc.
    if ($verbose){ $self->status_message("\n\n$rc_summary_cmd"); }
    mkdir($output_stats_dir);
    Genome::Sys->shellcmd(cmd => $rc_summary_cmd);
  }
 
  return 1;
}

1;

