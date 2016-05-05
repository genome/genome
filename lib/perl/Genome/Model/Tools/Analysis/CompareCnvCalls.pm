package Genome::Model::Tools::Analysis::CompareCnvCalls;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Analysis::CompareCnvCalls {
  is => 'Command::V2',
  has_input => [
    tp_bed => {
      is => 'FilesystemPath',
      doc => 'BED file containing true positive CNV calls',
    },
    eval_bed => {
      is => 'FilesystemPath',
      doc => 'BED file containg calls that need to be evaluated',
    },
    window_size => {
      is => 'Integer',
      doc => 'Size of the window to evaluate the CNV calls over',
      default => 500,
      is_optional => 1,
    },
    fai_file => {
      is => 'FilesystemPath',
      doc => '.fai index for the reference sequence, contains the length of the chromosomes',
      default => '/gscmnt/gc4096/info/model_data/2857786885/build102671028/all_sequences.fa.fai',
    },
    outdir => {
      is => 'FilesystemPath',
      doc => 'Directory to write results',
    },
    ROI_file => {
      is => 'FilesystemPath',
      doc => 'If the CNV calls were made on a specific region provide the bed file of this region.' . 
        ' The evaluations will be done on this region alone.',
      is_optional => 1,
    },
    sample => {
      is => 'Text',
      doc => 'Name of the sample',
      default => 'test',
      is_optional => 1,
    },
    test => {
      is => 'Boolean',
      doc => 'True for tests',
      default_value => 0,
      is_optional => 1,
    },
  ],
  doc => 'Compare CNV calls in two BED files.'
};

sub help_synopsis {
  return <<EOS
        gmt analysis compare-cnv-calls --outdir=/gscuser/gscuser1/tmp/
        --tp-bed=tp.bed --eval-bed=eval.bed
EOS
}

sub help_detail {
  return <<EOS
Compare CNV calls in two BED files. One of the files - tp-bed - is assumed to be the true-positive CNV calls. The other file - eval-bed - contains the calls that you would like to evaluate. The tool works by splitting the genome into fixed-width non-overlapping windows and classifying each window and classifying each window into either a TP or TN or NA based on tp-bed. A window is called a TP if it completely falls inside a CNV region in eval-bed AND completely falls inside a CNV region in the tp-bed file. A window is called a TN if it completely falls outside a CNV region in eval-bed AND it completely falls outside a CNV region in the tp-bed file. A window is called a FP if it completely falls inside a CNV region in eval-bed AND it completely falls outside a CNV region in the tp-bed file. A window is called a FN if it completely falls outside a CNV region in eval-bed AND it completely falls inside a CNV region in the tp-bed file. Any windows that have partial overlap with the CNV regions in eval-bed OR tp-bed are classified as NA. This approach has been well described in the cn.mops paper - "cn.MOPS: mixture of Poissons for discovering copy number variations in next-generation sequencing data with a low false discovery rate.
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

sub get_chr_sizes {
  my $self = shift;
  my $chr_size = shift;
  my $fai_file = $self->fai_file;
  open(my $fai_fh, $fai_file) 
    or die Genome::Sys->error_message("Unable to open $fai_file");
  while(<$fai_fh>) {
    my @splits = split(/\s+/, $_);
    my $chr = $splits[0];
    my $size = $splits[1];
    $chr_size->{$chr} = $size;
  }
  close($fai_fh);
}

sub copy_sort_bed {
  my $self = shift;
  my $bed = shift;
  my $out_bed_sorted = $self->outdir . "/" . basename($bed) . ".sorted";
  Genome::Sys->shellcmd(cmd => "joinx sort $bed -o $out_bed_sorted");
  return $out_bed_sorted;
}

#these files are assumed to be sorted
sub joinx_intersect {
  my $self = shift;
  my $bed_a = shift;
  my $bed_b = shift;
  my $bed_op = shift;
  my $bed_miss_a = shift;
  my $joinx_intersect = "joinx intersect $bed_a $bed_b -o $bed_op --miss-a $bed_miss_a";
  Genome::Sys->shellcmd(cmd=>$joinx_intersect);
}

sub write_ROC_metrics {
  my $self = shift;
  my $metrics_f = shift;
  my $cumul_metrics = shift;
  my $outdir = shift;
  open(my $METRICS_FH, ">", $metrics_f);
  print $METRICS_FH "sample\tTotal_P_windows\tTotal_N_windows\tTP_windows\tTN_windows\tFP_windows\tFN_windows\tTPR_Sensitivity".
    "\tTNR\tFPR\tPPV\n";
  foreach my $sample (keys %$cumul_metrics) {
    my ($TPR, $TNR, $FPR, $PPV);
    $TPR = eval ('$cumul_metrics->{$sample}{"TP_windows"} /
      ($cumul_metrics->{$sample}{"TP_windows"} + $cumul_metrics->{$sample}{"FN_windows"})');
    if($@) {
      $TPR = 0;
    }
    $TNR = eval ('$cumul_metrics->{$sample}{"TN_windows"} /
      ($cumul_metrics->{$sample}{"TN_windows"} + $cumul_metrics->{$sample}{"FP_windows"})');
    if($@) {
      $TNR = 0;
    }
    $FPR = eval ('$cumul_metrics->{$sample}{"FP_windows"} /
      ($cumul_metrics->{$sample}{"FP_windows"} + $cumul_metrics->{$sample}{"TP_windows"})');
    if($@) {
      $FPR = 0;
    }
    $PPV = eval ('$cumul_metrics->{$sample}{"TP_windows"} /
      ($cumul_metrics->{$sample}{"TP_windows"} + $cumul_metrics->{$sample}{"FP_windows"})');
    if($@) {
      $PPV = 0;
    }
    print $METRICS_FH $sample . "\t" .
    $cumul_metrics->{$sample}{"P_windows"} . "\t" .
    $cumul_metrics->{$sample}{"N_windows"} . "\t" .
    $cumul_metrics->{$sample}{"TP_windows"} . "\t" .
    $cumul_metrics->{$sample}{"TN_windows"} . "\t" .
    $cumul_metrics->{$sample}{"FP_windows"} . "\t" .
    $cumul_metrics->{$sample}{"FN_windows"} .  "\t" .
    $TPR . "\t" . $TNR . "\t" . $FPR . "\t" . $PPV .
    "\n";
  }
  close($METRICS_FH);
}

sub accumulate_ROC_metrics {
  my $self = shift;
  my $sample = shift;
  my $cumul_metrics = shift;
  my $outdir = shift;
  my $p_windows = $outdir . "/" . $sample . ".p.windows.bed";
  my $tp_windows = $outdir . "/" . $sample . ".tp.windows.bed";
  my $n_windows = $outdir . "/" . $sample . ".n.windows.bed";
  my $tn_windows = $outdir . "/" . $sample . ".tn.windows.bed";
  my $fn_windows = $outdir . "/" . $sample . ".fn.windows.bed";
  my $fp_windows = $outdir . "/" . $sample . ".fp.windows.bed";
  my $p_c = `wc -l < $p_windows`;
  my $n_c = `wc -l < $n_windows`;
  my $TP_windows = `wc -l < $tp_windows`;
  my $TN_windows = `wc -l < $tn_windows`;
  my $FP_windows = `wc -l < $fp_windows`;
  my $FN_windows = `wc -l < $fn_windows`;
  chomp ($p_c, $n_c, $TP_windows, $TN_windows, $FP_windows, $FN_windows);
  $cumul_metrics->{$sample}{"P_windows"} = $p_c;
  $cumul_metrics->{$sample}{"N_windows"} = $n_c;
  $cumul_metrics->{$sample}{"TP_windows"} = $TP_windows;
  $cumul_metrics->{$sample}{"TN_windows"} = $TN_windows;
  $cumul_metrics->{$sample}{"FP_windows"} = $FP_windows;
  $cumul_metrics->{$sample}{"FN_windows"} = $FN_windows;
} 

sub calculate_stats {
  my $self = shift;
  my $window_file = shift;
  my $sample = $self->sample;
  my $tp_bed_sorted = $self->copy_sort_bed($self->tp_bed);
  my $eval_bed_sorted = $self->copy_sort_bed($self->eval_bed);
  my $p_windows = $self->outdir . "/" . $sample . ".p.windows.bed";
  my $tp_windows = $self->outdir . "/" . $sample . ".tp.windows.bed";
  my $n_windows = $self->outdir . "/" . $sample . ".n.windows.bed";
  my $tn_windows = $self->outdir . "/" . $sample . ".tn.windows.bed";
  my $fn_windows = $self->outdir . "/" . $sample . ".fn.windows.bed";
  my $fp_windows = $self->outdir . "/" . $sample . ".fp.windows.bed";
  $self->joinx_intersect($window_file,$tp_bed_sorted, $p_windows, $n_windows);
  $self->joinx_intersect($p_windows,$eval_bed_sorted, $tp_windows, $fn_windows);
  $self->joinx_intersect($n_windows,$eval_bed_sorted, $fp_windows, $tn_windows);
}

sub create_window_file {
  my $self = shift;
  my $window_size = $self->window_size;
  if($window_size < 1) {
    die $self->error_message("Enter a valid window size");
  }
  my $window_file = $self->outdir . "/windows.$window_size.bed";
  if(-e $window_file) {
    $self->status_message("Using existing window file $window_file");
    return $window_file;
  }
  my %chr_size;
  $self->get_chr_sizes(\%chr_size);
  my @chrs;
  if($self->test) {
      @chrs = (1);
  }
  else {
    @chrs = (1..22);
    push(@chrs, "X");
    push(@chrs, "Y");
  }
  foreach my $chr (@chrs) {
    if($chr !~ /GL/ and $chr !~ /MT/) {
      my $size = $chr_size{$chr};
      my $awk_cmd = "awk 'BEGIN { pos = 1; while(pos < $size) " . 
          "{ print \"$chr\\t\"pos\"\\t\"pos+$window_size-1 >> \"$window_file\"; pos += $window_size;}}'";
      print $awk_cmd;
      Genome::Sys->shellcmd(cmd=>$awk_cmd);
    }
  }
  if($self->ROI_file) {
    my $sorted_ROI = $self->copy_sort_bed($self->ROI_file);
    my $intersected_file = $window_file . ".intersectedwithROI";
    my ($miss_file_fh, $miss_file) = Genome::Sys->create_temp_file();
    $self->joinx_intersect($window_file, $sorted_ROI, $intersected_file, $miss_file);
    $window_file = $intersected_file;
  }
  return $window_file;
}

sub execute {
  my $self = shift;
  my $window_file = $self->create_window_file();
  $self->calculate_stats($window_file);
  return 1;
}

1;
