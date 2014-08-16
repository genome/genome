package Genome::Model::ClinSeq::Command::Converge::PlotCnvs;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::PlotCnvs {
  is => 'Genome::Model::ClinSeq::Command::Converge::Base',
  has_input => [
  calculate_metrics => {
    is => 'Boolean',
    doc => 'Flag to calculate ROC metrics.',
    is_optional => 1,
    default => 0,
  },
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory to write results',
  },
  ],
  doc => 'Plot CNV calls integrated from WGS, Exome and Microarrays.',
};

sub help_synopsis {
  return <<EOS
genome model clin-seq converge plot-cnvs --builds='id in ["4b7539bb10cc4b9c97577cf11f4c79a2","cdca0edf526c4fe193d3054627a5871b"]' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge plot-cnvs --builds='model.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

genome model clin-seq converge plot-cnvs --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/ --calculate_metrics

genome model clin-seq converge plot-cnvs --builds='model.id in ["279f50e35d2b479ea3c32486eafd4ad4","7143119a93984056ae3f32c88c9ac2a1"],is_last_complete=1' --outdir=/tmp/snv_indel_report/
EOS
}

sub help_detail {
  return <<EOS
Plot CNVs from different technologies(microarray, exome etc), helps contrast and evaluate the calls.
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);

  unless (-e $self->outdir && -d $self->outdir) {
    push @errors, UR::Object::Tag->create(
      type => 'error',
      properties => ['outdir'],
      desc => "Outdir: " . $self->outdir . " not found or not a directory",
    );
  }
  return @errors;
}

sub copy_to_outdir {
  my $self = shift;
  my $file = shift;
  my $outfile = shift;
  my $outdir = $self->outdir;
  $outfile = $outdir . "/" . $outfile;
  unlink $outfile;
  Genome::Sys->copy_file($file, $outfile);
  return $outfile;
}

sub copy_wg_plots {
  my $self = shift;
  my $clinseq_build = shift;
  my $common_name = shift;
  my $microarray_wg_plot = $clinseq_build->microarray_cnv_wg_plot();
  if($microarray_wg_plot) {
    $microarray_wg_plot = $self->copy_to_outdir($microarray_wg_plot, $common_name . ".microarrays.cnvs.jpg");
  }
  my $wgs_wg_plot = $clinseq_build->wgs_cnv_wg_plot();
  if($wgs_wg_plot) {
    $wgs_wg_plot = $self->copy_to_outdir($wgs_wg_plot, $common_name . ".wgs.cnvs.jpg");
  }
  my $exome_wg_plot = $clinseq_build->exome_cnv_wg_plot();
  if($exome_wg_plot) {
    $exome_wg_plot = $self->copy_to_outdir($exome_wg_plot, $common_name . ".exome.cnvs.pdf");
  }
}

sub copy_wg_bed {
  my $self = shift;
  my $clinseq_build = shift;
  my $common_name = shift;
  my $microarray_cnvhmm_file = $clinseq_build->microarray_cnvhmm_file();
  if($microarray_cnvhmm_file) {
    $microarray_cnvhmm_file = $self->copy_to_outdir($microarray_cnvhmm_file, $common_name . ".microarrays.cnvs.txt");
  } else {
    $microarray_cnvhmm_file = "NA";
  }
  my $wgs_cnvhmm_file = $clinseq_build->wgs_cnvhmm_file();
  if($wgs_cnvhmm_file) {
    $wgs_cnvhmm_file = $self->copy_to_outdir($wgs_cnvhmm_file, $common_name . ".wgs.cnvs.txt");
  } else {
    $wgs_cnvhmm_file = "NA";
  }
  my $exome_cnvs_file = $clinseq_build->exome_cnvs_file();
  if($exome_cnvs_file) {
    $exome_cnvs_file = $self->copy_to_outdir($exome_cnvs_file, $common_name . ".exome.cnvs.txt");
  } else {
    $exome_cnvs_file = "NA";
  }
  return ($microarray_cnvhmm_file, $wgs_cnvhmm_file, $exome_cnvs_file); 
}

sub copy_files {
  my $self = shift;
  my $build = shift;
  my $common_name = shift;
  my $clinseq_build = Genome::Model::Build::ClinSeq->get(id => $build->id);
  $self->copy_wg_plots($clinseq_build, $common_name);
  return $self->copy_wg_bed($clinseq_build, $common_name);
}

sub get_cn {
  my $self = shift;
  my $ip_f = shift;
  my $op_f = shift;
  my $format = "awk '!/chr|CHR|GL/ { print \$1\"\\t\"\$2\"\\t\"\$3\"\\tNA\\t\"\$7-\$9+2 }' $ip_f > $op_f";
  Genome::Sys->shellcmd(cmd => $format);
  return;
}

sub get_bed {
  my $self = shift;
  my $ip_f = shift;
  my $op_f = shift;
  #insert a '0 0 1' line to deal with empty files.
  my $format = "awk '!/chr|CHR/ { print \$1\"\\t\"\$2\"\\t\"\$3 }' $ip_f > $op_f";
  Genome::Sys->shellcmd(cmd => $format);
  return;
}

sub get_exome_bed {
  my $self = shift;
  my $ip_f = shift;
  my $op_f = shift;
  my $format = "awk '!/chr|CHR/ { if(\$5 > 0.5) { print \$1\"\\t\"\$2\"\\t\"\$3 } }' $ip_f > $op_f";
  Genome::Sys->shellcmd(cmd => $format);
  return;
}

sub format_files {
  my $self = shift;
  my $microarray_file = shift;
  my $wgs_file = shift;
  my $exome_file = shift;
  my $microarray_cn = $$microarray_file . ".cn";
  my $wgs_cn = $$wgs_file . ".cn";
  my $exome_cn = $$exome_file . ".cn";
  my $microarray_bed = $$microarray_file . ".bed";
  my $wgs_bed = $$wgs_file . ".bed";
  my $exome_bed = $$exome_file . ".bed";
  $self->get_cn($$microarray_file, $microarray_cn);
  $self->get_cn($$wgs_file, $wgs_cn);
  $self->get_cn($$exome_file, $exome_cn);
  $self->get_bed($$microarray_file, $microarray_bed);
  $self->get_bed($$wgs_file, $wgs_bed);
  $self->get_exome_bed($$exome_file, $exome_bed);
  $self->joinxSortFile(\$microarray_bed);
  $self->joinxSortFile(\$wgs_bed);
  $self->joinxSortFile(\$exome_bed);
  $$microarray_file = $microarray_cn;
  $$exome_file = $exome_cn;
  $$wgs_file = $wgs_cn;
  return ($microarray_bed, $wgs_bed, $exome_bed);
}

sub create_combined_plot {
  my $self = shift;
  my $microarray_file = shift;
  my $wgs_file = shift;
  my $exome_file = shift;
  my $common_name = shift;
  my $plot = $self->outdir . $common_name . ".combinedCNV.pdf";
  my $Rscript = "Rscript " . __FILE__ . ".R";
  my $plot_cmd = $Rscript . " $microarray_file $wgs_file $exome_file $plot $common_name";
  Genome::Sys->shellcmd(cmd => $plot_cmd);
}

sub calculate_ROC_metrics {
  my $self = shift;
  my $tp_bed = shift;
  my $eval_bed = shift;
  my $common_name = shift;
  my $outdir = $self->outdir;
  my $c1 = Genome::Model::Tools::Analysis::CompareCnvCalls->create(
    tp_bed => $tp_bed, eval_bed => $eval_bed, outdir => $outdir,
    sample => $common_name);
  $c1->execute();
}

sub joinxSortFile {
  my ($self, $ip_file) = @_;
  my $sorted_op_file = $$ip_file . ".sorted";
  my $joinx_sort_cmd = Genome::Model::Tools::Joinx::Sort->create(output_file=>$sorted_op_file, input_files=>[$ip_file]);
  $joinx_sort_cmd->execute();
  unlink $$ip_file;
  $$ip_file = $sorted_op_file;
}

sub plot_wgs_exome_microarray_cnvs() {
  my $self = shift;
  my @builds = $self->builds;
  my (%cumulative_metrics1, %cumulative_metrics2, %cumulative_metrics3);
  foreach my $build (@builds) {
    my $microarray_file;
    my $wgs_file;
    my $exome_file;
    my $common_name = $build->common_name;
    my $subject_common_name = $build->model->exome_model->tumor_model->subject->common_name;
    $subject_common_name =~ s/[\s,-]/_/g;
    $common_name = $common_name . "_" . $subject_common_name;
    $common_name =~ s/_tumor//;
    $common_name =~ s/_relapse2//;
    $common_name =~ s/_pre_treatment_met/_1/;
    $common_name =~ s/_recurrence_met/_2/;
    ($microarray_file, $wgs_file, $exome_file) = $self->copy_files($build, $common_name);
    if ($microarray_file eq "NA" or $wgs_file eq "NA" or $exome_file eq "NA") {
      $self->status_message("Skipping $common_name, this sample does not have all three CNV files.");
      next;
    }
    my ($microarray_bed, $wgs_bed, $exome_bed) = $self->format_files(
      \$microarray_file, \$wgs_file, \$exome_file);
    if($self->calculate_metrics) {
      $self->calculate_ROC_metrics($wgs_bed, $exome_bed, $common_name);
      $self->calculate_ROC_metrics($wgs_bed, $microarray_bed, $common_name);
      $self->calculate_ROC_metrics($exome_bed, $microarray_bed, $common_name);
      Genome::Model::Tools::Analysis::CompareCnvCalls->accumulate_ROC_metrics(
        $common_name, \%cumulative_metrics1, $self->outdir);
      Genome::Model::Tools::Analysis::CompareCnvCalls->accumulate_ROC_metrics(
        $common_name, \%cumulative_metrics2, $self->outdir);
      Genome::Model::Tools::Analysis::CompareCnvCalls->accumulate_ROC_metrics(
        $common_name, \%cumulative_metrics3, $self->outdir);
    }
    $self->create_combined_plot($microarray_file, $wgs_file, $exome_file,
      $common_name);
  }
  if($self->calculate_metrics) {
    my $metrics_f1 = $self->outdir . "/" . "wgs_exome.ROC.metrics.txt";
    Genome::Model::Tools::Analysis::CompareCnvCalls->write_ROC_metrics(
      $metrics_f1, \%cumulative_metrics1, $self->outdir);
    my $metrics_f2 = $self->outdir . "/" . "wgs_ma.ROC.metrics.txt";
    Genome::Model::Tools::Analysis::CompareCnvCalls->write_ROC_metrics(
      $metrics_f2, \%cumulative_metrics2, $self->outdir);
    my $metrics_f3 = $self->outdir . "/" . "exome_ma.ROC.metrics.txt";
    Genome::Model::Tools::Analysis::CompareCnvCalls->write_ROC_metrics(
      $metrics_f3, \%cumulative_metrics3, $self->outdir);
  }
}

sub execute {
  my $self = shift;
  $self->plot_wgs_exome_microarray_cnvs();
  return 1;
}

1;
