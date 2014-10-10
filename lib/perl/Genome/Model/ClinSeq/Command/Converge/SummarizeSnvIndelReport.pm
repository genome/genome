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
  },
  outdir => {
    is => 'FilesystemPath',
    doc => 'Directory to write results',
  },
  min_mq => {
    is => 'Integer',
    doc => 'minimum mapping quality of reads to be considered',
    default => '30',
  },
  min_bq => {
    is => 'Integer',
    doc => 'minimum base quality of bases in reads to be considered',
    default => '20',
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
  my $build = $self->clinseq_build;
  my $bq = $self->min_bq;
  my $mq = $self->min_mq;
  my $snv_indel_report = $build->snv_indel_report_clean_filtered_file($bq, $mq);
  my (%snv_stats, %indel_stats);
  $self->parse_snv_indel_report($snv_indel_report, \%snv_stats, \%indel_stats);
  $self->write_stats(\%snv_stats, \%indel_stats);
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
    foreach my $caller (@callers) {
      foreach my $data_type (@data_types) {
        if($stats->{$data_type}->{caller}->{$caller}) {
          $stats->{$data_type}->{caller}->{$caller} += 1;
        } else {
          $stats->{$data_type}->{caller}->{$caller} = 1;
        }
      }
    }
  }
}

sub write_stats {
  my $self = shift;
  my $snv_stats = shift;
  my $indel_stats = shift;
  my $stats_file = $self->outdir . "/Stats.tsv";
  open my $STATS, ">$stats_file";
  print $STATS "Question\tAnswer\tData_Type\tAnalysis_Type\tStatistic_Type\tExtra_Description\n";
  foreach my $data_source ("wgs", "exome") {
    $self->write_snv_stats($snv_stats, $data_source, $STATS);
    $self->write_indel_stats($indel_stats, $data_source, $STATS);
  }
  close $STATS;
}

sub write_snv_stats {
  my $self = shift;
  my $snv_stats = shift;
  my $data_source = shift;
  my $STATS = shift;
  foreach my $caller ("strelka", "sniper", "samtools", "varscan", "mutect") {
    unless ($snv_stats->{$data_source}->{caller}->{$caller}) {
      $snv_stats->{$data_source}->{caller}->{$caller} = 0;
    }
  }
  print $STATS "Number of Strelka SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"strelka"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Strelka\n";
  print $STATS "Number of Sniper SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"sniper"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Sniper\n";
  print $STATS "Number of VarScan SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"varscan"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by VarScan\n";
  print $STATS "Number of SamTools SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"samtools"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by SamTools\n";
  print $STATS "Number of Mutect SNV calls\t" . $snv_stats->{$data_source}->{caller}->{"mutect"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of SNVs called by Mutect\n";
}

sub write_indel_stats {
  my $self = shift;
  my $indel_stats = shift;
  my $data_source = shift;
  my $STATS = shift;
  foreach my $caller ("strelka", "gatk", "pindel", "varscan") {
    unless ($indel_stats->{$data_source}->{caller}->{$caller}) {
      $indel_stats->{$data_source}->{caller}->{$caller} = 0;
    }
  }
  print $STATS "Number of Strelka Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"strelka"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Strelka\n";
  print $STATS "Number of GATK Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"gatk"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of Indels called by GATK\n";
  print $STATS "Number of Pindel Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"pindel"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Pindel\n";
  print $STATS "Number of VarScan Indel calls\t" . $indel_stats->{$data_source}->{caller}->{"varscan"}.
  "\t" . $data_source . "\tClinseq Build Summary\tCount\tNumber of Indels called by Varscan\n";
}

1;
