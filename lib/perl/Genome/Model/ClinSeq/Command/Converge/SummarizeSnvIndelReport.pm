package Genome::Model::ClinSeq::Command::Converge::SummarizeSnvIndelReport;

use strict;
use warnings;
use Genome;

class Genome::Model::ClinSeq::Command::Converge::SummarizeSnvIndelReport {
  is => 'Command::V2',
  has_input => [
  clinseq_build => {
    is => 'Genome::Model::Build::ClinSeq',
    doc => 'ClinSeq build to summarize snv_indel_report for.',
    is_optional => 1,
  },
  filtered_report => {
    is => 'FilesystemPath',
    doc => 'Filtered SnvIndel report. If not specified, this file is obtained from the snv_indel report from the clinseq build.',
    is_optional => 1,
  },
  unfiltered_report => {
    is => 'FilesystemPath',
    doc => 'Unfiltered SnvIndel report. If not specified, this file is obtained from the snv_indel report from the clinseq build.',
    is_optional => 1,
  },
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory to write results',
  },
  min_mq => {
    is => 'Integer',
    doc => 'minimum mapping quality of reads to be considered',
  },
  min_bq => {
    is => 'Integer',
    doc => 'minimum base quality of bases in reads to be considered',
  },
  ],
  doc => 'Summarize SNVIndel Report in ClinSeq.',
};

sub help_synopsis {
  return <<EOS
        genome model clin-seq generate-sciclone-plots \\
        --outdir=/gscuser/gscuser1/tmp/ \\
        --clinseq-build='a4abcd1313eb4376b59e68a9dd9d5ad2'
EOS
}

sub help_detail {
  return <<EOS
Summarize the SNV-Indel Report in ClinSeq builds.
EOS
}

sub execute {
  my $self = shift;
  my $reports = $self->get_sireports();
  my $stats_file = $self->outdir . "/Stats.tsv";
  if(-e $stats_file) {
    unlink $stats_file;
  }
  foreach my $filter(keys %$reports) {
    my $report = $reports->{$filter};
    my (%snv_stats, %indel_stats);
    $self->parse_snv_indel_report($report, \%snv_stats, \%indel_stats);
    $self->write_stats(\%snv_stats, \%indel_stats, $stats_file, $filter);
  }
}

sub get_sireports {
  my $self = shift;
  my $reports;
  my $bq = $self->min_bq;
  my $mq = $self->min_mq;
  my ($snv_indel_report1, $snv_indel_report2);
  if($self->clinseq_build) {
    my $build = $self->clinseq_build;
    $snv_indel_report1 = $build->snv_indel_report_clean_filtered_file($bq, $mq);
    $snv_indel_report2 = $build->snv_indel_report_clean_unfiltered_file($bq, $mq);
  } elsif ($self->filtered_report and $self->unfiltered_report) {
    $snv_indel_report1 = $self->filtered_report;
    $snv_indel_report2 = $self->unfiltered_report;
  } else {
    die $self->error_message("SnvIndel reports or a clinseq build not specified.");
  }
  $reports->{"filtered"} = $snv_indel_report1;
  $reports->{"unfiltered"} = $snv_indel_report2;
  return $reports;
}

sub parse_snv_indel_report {
  my $self = shift;
  my $snv_indel_report = shift;
  my $snv_stats = shift;
  my $indel_stats = shift;
  my $stats;
  my $reader = Genome::Utility::IO::SeparatedValueReader->create(
    separator => "\t",
    input => $snv_indel_report,
  );
  while (my $data = $reader->next) {
    if ($data->{chromosome_name} =~ /MT|GL/) {
      next;
    }
    if ($data->{type} =~ /SNP/) {
      $stats = $snv_stats;
    } elsif ($data->{type} =~ /INS|DEL/) {
      $stats = $indel_stats;
    }
    my @callers = split(",", $data->{variant_source_callers});
    my @data_types = split(",", $data->{data_type});
    my @filters = split(",", $data->{filtered});
    foreach my $caller (@callers) {
      foreach my $data_type (@data_types) {
        if($stats->{$data_type}->{caller}->{$caller}) {
          $stats->{$data_type}->{caller}->{$caller} += 1;
        } else {
          $stats->{$data_type}->{caller}->{$caller} = 1;
        }
      }
    }
    foreach my $filter (@filters) {
      foreach my $data_type (@data_types) {
        if(defined $stats->{$data_type}->{filters}->{$filter}) {
          $stats->{$data_type}->{filters}->{$filter} += 1;
        } else {
          $stats->{$data_type}->{filters}->{$filter} = 1;
        }
      }
    }
  }
}

sub write_stats {
  my $self = shift;
  my $snv_stats = shift;
  my $indel_stats = shift;
  my $stats_file = shift;
  my $filter = shift;
  open my $STATS, ">>$stats_file";
  print $STATS "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";
  foreach my $data_source ("wgs", "exome") {
    $self->write_snv_stats($snv_stats, $data_source, $STATS, $filter);
    $self->write_indel_stats($indel_stats, $data_source, $STATS, $filter);
  }
  close $STATS;
}

sub write_snv_stats {
  my $self = shift;
  my $snv_stats = shift;
  my $data_source = shift;
  my $STATS = shift;
  my $filter = shift;
  foreach my $caller ("strelka", "sniper", "samtools", "varscan", "mutect") {
    unless ($snv_stats->{$data_source}->{caller}->{$caller}) {
      $snv_stats->{$data_source}->{caller}->{$caller} = 0;
    }
  }
  my $source = $data_source . "," . $filter;
  print $STATS "Number of Strelka SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"strelka"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Strelka\n";
  print $STATS "Number of Sniper SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"sniper"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Sniper\n";
  print $STATS "Number of VarScan SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"varscan"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by VarScan\n";
  print $STATS "Number of SamTools SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"samtools"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by SamTools\n";
  print $STATS "Number of Mutect SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"mutect"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Mutect\n";
  my $filters = $snv_stats->{$data_source}->{filters};
  print $STATS "Number of passing SNV calls\t" . eval{$snv_stats->{$data_source}->{filters}->{0}||0} .
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs passing all filters.\n";
  foreach my $filter(keys %$filters) {
    if($filter ne "0") {
      print $STATS "Number of SNV calls filtered out by $filter\t" . $snv_stats->{$data_source}->{filters}->{$filter} .
        "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of SNVs filtered by $filter\n";
    }
  }
}

sub write_indel_stats {
  my $self = shift;
  my $indel_stats = shift;
  my $data_source = shift;
  my $STATS = shift;
  my $filter = shift;
  foreach my $caller ("strelka", "gatk", "pindel", "varscan") {
    unless ($indel_stats->{$data_source}->{caller}->{$caller}) {
      $indel_stats->{$data_source}->{caller}->{$caller} = 0;
    }
  }
  my $source = $data_source . "," . $filter;
  print $STATS "Number of Strelka Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"strelka"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Strelka\n";
  print $STATS "Number of GATK Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"gatk"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of Indels called by GATK\n";
  print $STATS "Number of Pindel Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"pindel"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Pindel\n";
  print $STATS "Number of VarScan Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"varscan"}.
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Varscan\n";
  my $filters = $indel_stats->{$data_source}->{filters};
  print $STATS "Number of passing INDEL calls\t" . eval{$indel_stats->{$data_source}->{filters}->{0}||0} .
  "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of INDELs passing all filters.\n";
  foreach my $filter(keys %$filters) {
    if($filter ne "0") {
      print $STATS "Number of INDELs calls filtered out by $filter\t" . $indel_stats->{$data_source}->{filters}->{$filter} .
      "\t" . $source . "\tClinseq Build Summary\tCount\tNumber of INDELs filtered by $filter\n";
    }
  }
}

1;
