package Genome::Model::Tools::Analysis::CompareCnvCalls;

use strict;
use warnings;
use Genome;

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
      default => '/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa.fai',
    },
    outdir => {
      is => 'FilesystemPath',
      doc => 'Directory to write results',
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
        gmt copy-number analysis copy-number compare-cnv-calls --outdir=/gscuser/gscuser1/tmp/ --tp-bed=tp.bed --eval-bed=eval.bed
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
  my $out_bed = basename($bed);
  my $out_bed_sorted = basename($bed) . ".sorted";
  Genome::Sys->shellcmd(cmd => "cp $bed $out_bed");
  Genome::Sys->shellcmd(cmd => "joinx sort $out_bed -o $out_bed_sorted");
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

sub calculate_stats {
  my $self = shift;
  my $window_file = shift;
  my $tp_bed_sorted = $self->copy_sort_bed($self->tp_bed);
  my $eval_bed_sorted = $self->copy_sort_bed($self->eval_bed);
  my $window_file_sorted = $self->copy_sort_bed($self->window_file);
  my $p_windows = $self->outdir . "p.windows.bed";
  my $tp_windows = $self->outdir . "tp.windows.bed";
  my $n_windows = $self->outdir . "n.windows.bed";
  my $tn_windows = $self->outdir . "tn.windows.bed";
  my $fn_windows = $self->outdir . "fn.windows.bed";
  my $fp_windows = $self->outdir . "fp.windows.bed";
  $self->joinx_intersect($window_file_sorted, $tp_bed_sorted, $p_windows, $n_windows);
  $self->joinx_intersect($p_windows, $eval_bed_sorted, $tp_windows, $fp_windows);
  $self->joinx_intersect($n_windows, $eval_bed_sorted, $fp_windows, $tn_windows);
}

sub calculate_true_negative {
  my $self = shift;
}

sub calculate_false_negative {
  my $self = shift;
}

sub create_window_file {
  my $self = shift;
  my $window_size = $self->window_size;
  if($window_size < 1) {
    die $self->error_message("Enter a valid window size");
  }
  my $window_file = $self->outdir . "/windows.$window_size.bed";
  Genome::Sys->shellcmd(cmd => "rm -f $window_file");
  my %chr_size;
  $self->get_chr_sizes(\%chr_size);
  my @chrs = (1..22);
  push(@chrs, "X");
  push(@chrs, "Y");
  foreach my $chr (@chrs) {
    print $chr;
    if($chr !~ /GL/ and $chr !~ /MT/) {
      my $size = $chr_size{$chr};
      my $awk_cmd = "awk 'BEGIN { pos = 1; while(pos < $size) " . 
          "{ print \"$chr\\t\"pos\"\\t\"pos+$window_size-1 >> \"$window_file\"; pos += $window_size;}}'";
      print $awk_cmd;
      Genome::Sys->shellcmd(cmd=>$awk_cmd);
    }
  }
  return $window_file;
}

sub execute {
  my $self = shift;
  my $window_file = $self->create_window_file();
  $self->calculate_stats($window_file);
  print "window file is $window_file";
  return 1;
}

1;
