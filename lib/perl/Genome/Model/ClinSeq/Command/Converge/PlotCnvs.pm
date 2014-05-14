package Genome::Model::ClinSeq::Command::Converge::PlotCnvs;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::PlotCnvs {
  is => 'Genome::Model::ClinSeq::Command::Converge::Base',
  has_input => [
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

genome model clin-seq converge plot-cnvs --builds='model_groups.id=9d0fcdca2b5d4f4385b83d2f75addac4,is_last_complete=1' --outdir=/tmp/snv_indel_report/

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

sub get_models() {
  my $self = shift;
  my $clinseq_mg = $self->clinseq_mg;  
  my @models = $clinseq_mg->models;
  return @models;
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
  my $microarray_cnv_wg_plot = $clinseq_build->microarray_cnv_wg_plot();
  if($microarray_cnv_wg_plot) {
    $microarray_cnv_wg_plot = $self->copy_to_outdir($microarray_cnv_wg_plot, $common_name . ".microarrays.cnvs.jpg");
  }
  my $wgs_cnv_wg_plot = $clinseq_build->wgs_cnv_wg_plot();
  if($wgs_cnv_wg_plot) {
    $wgs_cnv_wg_plot = $self->copy_to_outdir($wgs_cnv_wg_plot, $common_name . ".wgs.cnvs.jpg");
  }
  my $exome_cnv_wg_plot = $clinseq_build->exome_cnv_wg_plot();
  if($exome_cnv_wg_plot) {
    $exome_cnv_wg_plot = $self->copy_to_outdir($exome_cnv_wg_plot, $common_name . ".exome.cnvs.pdf");
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

sub copy_cnv_files {
  my $self = shift;
  my $build = shift;
  my $common_name = shift;
  my $clinseq_build = Genome::Model::Build::ClinSeq->get(id => $build->id);
  $self->copy_wg_plots($clinseq_build, $common_name);
  return $self->copy_wg_bed($clinseq_build, $common_name);
}

sub format_cnv_files {
  my $self = shift;
  my $microarray_cnv_file = shift;
  my $wgs_cnv_file = shift;
  my $exome_cnv_file = shift;
  my $microarray_cnv_file_f = $microarray_cnv_file . ".f";
  my $wgs_cnv_file_f = $wgs_cnv_file . ".f";
  my $exome_cnv_file_f = $exome_cnv_file . ".f";
  my $ma_format_cmd = "awk '!/chr|CHR|GL/ { print \$1\"\\t\"\$2\"\\t\"\$3\"\\tNA\\t\"\$7-\$9+2 }' $microarray_cnv_file > $microarray_cnv_file_f";
  my $wgs_format_cmd = "awk '!/chr|CHR/ { print \$1\"\\t\"\$2\"\\t\"\$3\"\\tNA\\t\"\$7-\$9+2 }' $wgs_cnv_file > $wgs_cnv_file_f";
  my $exome_format_cmd = "awk '!/chr|CHR/ { print \$1\"\\t\"\$2\"\\t\"\$3\"\\tNA\\t\"2^\$5*2 }' $exome_cnv_file > $exome_cnv_file_f";
  Genome::Sys->shellcmd(cmd => $ma_format_cmd);
  Genome::Sys->shellcmd(cmd => $wgs_format_cmd);
  Genome::Sys->shellcmd(cmd => $exome_format_cmd);
  return ($microarray_cnv_file_f, $wgs_cnv_file_f, $exome_cnv_file_f);
}

sub create_combined_plot {
  my $self = shift;
  my $microarray_cnv_file = shift;
  my $wgs_cnv_file = shift;
  my $exome_cnv_file = shift;
  my $common_name = shift;
  my $plot = $self->outdir . $common_name . ".combinedCNV.pdf";
  my $Rscript = "Rscript " . __FILE__ . ".R";
  my $plot_cmd = $Rscript . " $microarray_cnv_file $wgs_cnv_file $exome_cnv_file $plot $common_name";
  Genome::Sys->shellcmd(cmd => $plot_cmd);
}

sub plot_wgs_exome_microarray_cnvs() {
  my $self = shift;
  my @builds = $self->builds;
  foreach my $build (@builds) {
      my $microarray_cnv_file;
      my $wgs_cnv_file;
      my $exome_cnv_file;
      my $common_name = $build->common_name;
      ($microarray_cnv_file, $wgs_cnv_file, $exome_cnv_file) = $self->copy_cnv_files($build, $common_name);
      if ($microarray_cnv_file eq "NA" or $wgs_cnv_file eq "NA" or $exome_cnv_file eq "NA") {
        $self->status_message("Skipping $common_name, this sample does not have all three CNV files.");
        next;
      }
      ($microarray_cnv_file, $wgs_cnv_file, $exome_cnv_file) = $self->format_cnv_files($microarray_cnv_file, $wgs_cnv_file, $exome_cnv_file);
      $self->create_combined_plot($microarray_cnv_file, $wgs_cnv_file, $exome_cnv_file, $common_name);
  }
}

sub execute {
  my $self = shift;
  $self->plot_wgs_exome_microarray_cnvs();
  return 1;
}

1;
