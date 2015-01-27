package Genome::Model::ClinSeq::Command::IdentifyLoh;

use strict;
use warnings;
use File::Spec;
use Genome;

class Genome::Model::ClinSeq::Command::IdentifyLoh {
  is => ['Command::V2',
  'Genome::Model::ClinSeq::Util',
  ],
  has_input => [
  clinseq_build => {
    is => 'Genome::Model::Build::ClinSeq',
    doc => 'ClinSeq build to identify LOH regions in.
      [Either this or a somatic variation build is required.]',
    is_optional => 1,
  },
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory where output files will be written.',
  },
  somvar_build => {
    is => 'Genome::Model::Build::ClinSeq',
    doc => 'SomVar build to identify LOH regions in.
      [Either this or a clinseq build is required.]',
    is_optional => 1,
  },
  minprobes=> {
    is => 'Number',
    is_optional => 1,
    doc => 'Minimum number of probes to call a loh segment.',
    default => 10,
  },
  segment_pc_cutoff => {
    is => 'Number',
    is_optional => 1,
    doc => 'Minimum percent of LOH in a segment. n_loh_probes/n_probes * 100',
    default => 0.95,
  },
  test => {
    is => 'Boolean',
    doc => 'set for test-cases',
    is_optional => 1,
    default => 0,
  },
  ],
  doc => 'Identify regions of LOH using clinseq tumor/normal pairs.',
};

sub help_synopsis {
  return <<EOS
genome model clin-seq identify-loh\\
  --outdir=/gscuser/gscuser1/tmp/ \\
  --clinseq-build='a4abcd1313eb4376b59e68a9dd9d5ad2'
genome model clin-seq identify-loh\\
  --outdir=/gscuser/gscuser1/tmp/ \\
  --somvar-build='a4abcd1313eb4376b59e68a9dd9d5ad2'
EOS
}

sub help_detail {
  return <<EOS
Use the results from Varscan in the somatic-variation builds and identify
regions of LOH. Specifically, look at heterozygous single nucleotide variants
in the normal sample and compare with stretches of homozygosity in the tumor.
We will attempt to use the results from the WGS-somatic-variation build input, if
these are unavailable the tool looks for WEx-somatic-variation build.

Alternatively, instead of a clin-seq model you can supply a
somatic-variation model as an input.

Thanks to Chris Miller for working out the workflow.
EOS
}

sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__(@_);
  unless (-e $self->outdir && -d $self->outdir) {
    push @errors, UR::Object::Tag->create(
      type => 'error',
      properties => ['outdir'],
      desc => "Outdir: " . $self->outdir .
      " not found or not a directory",
    );
  }
  return @errors;
}

sub resolve_somvar {
  my $self = shift;
  my $somvar_build = $self->somvar_build;
  unless($somvar_build) {
    my $clinseq_build = $self->clinseq_build;
    unless($clinseq_build) {
      die $self->error_message("Please pass somvar or clinseq build as input.");
    }
    $somvar_build = $self->get_best_somvar_build($clinseq_build);
  }
  return $somvar_build;
}

sub split_snvs {
  my $self = shift;
  my $varscan_op = shift;
  my $snv_prefix = shift;
  my $awk_cmd;
  #use chr1 variants only for chr1
  if($self->test) {
    $awk_cmd = "awk '\$1==1 { print > \"$snv_prefix.\"\$1\".snp\" }' $varscan_op";
  } else {
    $awk_cmd = "awk '{ print > \"$snv_prefix.\"\$1\".snp\" }' $varscan_op";
  }
  Genome::Sys->shellcmd(cmd => $awk_cmd);
  return;
}

sub get_varscan_snvs {
  my $self = shift;
  my $somvar_build = shift;
  my $snv_prefix = shift;
  my $dd = $somvar_build->data_directory;
  my @varscan_snvs =
    glob(File::Spec->join($dd, "variants/snv/varscan-somatic-*/snvs.hq.unfiltered"));
  unless(scalar @varscan_snvs) {
    die $self->error_message("Unable to find varscan SNVs in $dd");
  }
  $self->status_message("varscan snvs file found: " .
    join(@varscan_snvs, "\t"));
  $self->split_snvs($varscan_snvs[0], $snv_prefix);
}

sub get_tumor_normal_bam {
  my $self = shift;
  my $somvar_build = shift;
  my $tumor_bam = $somvar_build->tumor_model->last_succeeded_build->
    whole_rmdup_bam_file;
  my $normal_bam = $somvar_build->normal_model->last_succeeded_build->
    whole_rmdup_bam_file;
  unless (-e $tumor_bam and -e $normal_bam) {
    die $self->error_message("tumor_bam $tumor_bam or normal_bam " .
      "$normal_bam does not exist.");
  }
  return ($tumor_bam, $normal_bam);
}

sub reference_build {
  my $self = shift;
  my $somvar_build = shift;
  my $normal_refalign_build = $somvar_build->normal_build if ($somvar_build);
  my $tumor_refalign_build = $somvar_build->tumor_build if ($somvar_build);
  my @input_builds = ($normal_refalign_build, $tumor_refalign_build, $somvar_build);
  for my $build (@input_builds){
    next unless $build;
    my $m = $build->model;
    if ($m->can("reference_sequence_build")){
      return $m->reference_sequence_build;
    }
  }
}

sub filter_snvs {
  my $self = shift;
  my $somvar_build = shift;
  my $snv_prefix = shift;
  my ($tumor_bam, $normal_bam) =
    $self->get_tumor_normal_bam($somvar_build);
  my $refbuild = $self->reference_build($somvar_build);
  my $ref_fa = $refbuild->full_consensus_path('fa');
  my $filter = Genome::Model::Tools::Varscan::SomaticParallelFilter->
    create(
      tumor_bam => $tumor_bam,
      normal_bam => $normal_bam,
      output => $snv_prefix,
      filter_germline => 1,
      reference => $ref_fa,
    );
  $filter->execute();
  Genome::Sys->shellcmd(cmd => "/gscuser/cmiller/usr/bin//bwait -n varscan -s 60");
}

sub combine_sort_snvs {
  my $self = shift;
  my $snv_prefix = shift;
  my $snv_combined = $snv_prefix . ".combined";
  my $snv_combined_sorted = $snv_prefix . ".combined.sorted";
  my $cat_cmd = "cat $snv_prefix.*.snp.*.LOH.hc.fpfilter " .
    "$snv_prefix.*.snp.*.Germline.hc.fpfilter " .
    "> $snv_combined";
  Genome::Sys->shellcmd(cmd => $cat_cmd);
  my $sort = Genome::Model::Tools::Capture::SortByChrPos->
    create(
      input_file => $snv_combined,
      output_file => $snv_combined_sorted,
    );
  $sort->execute;
  return $snv_combined_sorted;
}

sub segment_loh {
  my $self = shift;
  my $snvs_combined_sorted = shift;
  my $loh_basename = shift;
  my $segmenter = Genome::Model::Tools::Varscan::LohSegments->
    create(
      variant_file => $snvs_combined_sorted,
      output_basename => $loh_basename,
    );
  $segmenter->execute();
}

sub filter_loh {
  my $self = shift;
  my $loh_basename = shift;
  my $segment_pc_cutoff = $self->segment_pc_cutoff;
  my $minprobes = $self->minprobes;
  my $loh_segments =
    $loh_basename . ".segments.cbs";
  my $loh_segments_filtered = $loh_segments . ".filtered";
  unless(-e $loh_segments) {
    die $self->error_message("Unable to find loh file $loh_segments");
  }
  my $awk_filter =
    "awk '\$4 > $segment_pc_cutoff && \$3 >= $minprobes' $loh_segments " .
      "> $loh_segments_filtered";
  Genome::Sys->shellcmd(cmd => $awk_filter);
}

sub cleanup {
  my $self = shift;
  my $unwanted;
  foreach my $suffix ("*.Somatic*", "*.readcounts", "*.hc", "*.LOH",
    "*.lc", "*.Germline", "*.removed", "*.formatted",
    "*.hc.err", "*.other") {
    my $unwanted = File::Spec->join($self->outdir, $suffix);
    if(-e $unwanted) {
      unlink $unwanted;
    }
  }
}

sub execute {
  my $self = shift;
  my $somvar_build = $self->resolve_somvar;
  my $snv_prefix = File::Spec->join(
      $self->outdir,
      "varscan.snps");
  $self->get_varscan_snvs($somvar_build, $snv_prefix);
  $self->filter_snvs($somvar_build, $snv_prefix);
  my $combined_sorted = $self->combine_sort_snvs($snv_prefix);
  my $loh_basename = File::Spec->join(
    $self->outdir,
    "varscan.output.loh");
  my $loh_segments = $self->segment_loh($combined_sorted, $loh_basename);
  $self->filter_loh($loh_basename);
  $self->cleanup;
  return 1;
}

1;

