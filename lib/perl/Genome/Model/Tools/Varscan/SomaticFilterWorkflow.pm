package Genome::Model::Tools::Varscan::SomaticFilterWorkflow;

use strict;
use warnings;

use Genome;
use Workflow::Simple;
use File::Basename;

BEGIN {
    $ENV{WF_USE_FLOW} = 1;
}

class Genome::Model::Tools::Varscan::SomaticFilterWorkflow {
    is => 'Genome::Model::Tools::Varscan',
    has_input => [
        normal_bam => {
            is => 'Text',
            doc => "Path to Normal BAM file",
        },
        tumor_bam => {
            is => 'Text',
            doc => "Path to Tumor BAM file",
        },
        prefix => {
            is => 'Text',
            doc => "Input VarScan SNV prefix. The SNV files are found by looking for \$prefix.\$chr.snp",
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs",
            example_values => [(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')],
        },
        chromosome => {
            is => 'Text',
            doc => "Specify a single chromosome (optional)",
            is_optional => 1,
        },
        filter_loh => {
            is => 'Text',
            doc => "Flag to apply filter to LOH-HC calls using normal BAM",
            is_optional => 1,
            default => 1,
        },
        filter_germline => {
            is => 'Text',
            doc => "Flag to apply filter to Germline-HC calls using tumor BAM",
            is_optional => 1,
            default => 0,
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "If set to 1, skip execution if output files exist",
            is_optional => 1,
        },
        varscan_params => {
            is => 'Text',
            doc => "Parameters to pass to VarScan [--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1]" ,
            is_optional => 1,
        },
        outdir => {
            is => 'FilesystemPath',
            doc => 'Directory where output files will be written',
        },
    ],
};

sub help_brief {
    "This filter is identical to SomaticParallelFilter. Uses a workflow instead of manual 'bsub'. Filters SNVs only."
}

sub help_synopsis {
    return <<EOS
Filter Varscan somatic SNVs
gmt varscan somatic-parallel-filter --normal-bam nbam --tumor-bam  tbam --output varscan.snps
EOS
}

sub get_variant_files {
  my $self = shift;
  my $variant_files = shift;
  my $prefix = $self->prefix;
  my $reference = $self->reference;
  my $index_file = "$reference.fai";
  unless(-e $index_file) {
    die $self->error_message("Index file for reference ($index_file) not found!\n");
  }
  my $input = new FileHandle ($index_file);
  while (<$input>) {
    chomp;
    my $line = $_;
    my ($chrom) = split(/\t/, $line);
    next if($chrom =~ 'NT_' || $chrom =~ /GL/ || $chrom =~ /MT/);
    my $variant_f = $prefix . ".$chrom.snp";
	if(-e $variant_f) {
      push @$variant_files, $variant_f;
    }
  }
}

sub map_workflow_inputs {
  my $self = shift;
  my $input = shift;
  my @variant_files;
  $self->get_variant_files(\@variant_files);
  foreach my $vf(@variant_files) {
    push @{$input->{variants_file}}, $vf;
    my $formatted_output = File::Spec->join(
      $self->outdir,
      basename($vf) . ".formatted",
    );
    my $processed_somatic = File::Spec->join(
      $self->outdir,
      basename($vf) . ".formatted.Somatic.hc",
    );
    my $processed_germline = File::Spec->join(
      $self->outdir,
      basename($vf) . ".formatted.Germline.hc",
    );
    my $processed_loh = File::Spec->join(
      $self->outdir,
      basename($vf) . ".formatted.LOH.hc",
    );
    my $unfiltered_output = File::Spec->join(
      $self->outdir,
      basename($vf) . ".formatted.unfiltered",
    );
    push @{$input->{format_output}}, $formatted_output;
    push @{$input->{unfiltered}}, $unfiltered_output;
    push @{$input->{processed_somatic}}, $processed_somatic;
    push @{$input->{processed_germline}}, $processed_germline;
    push @{$input->{processed_loh}}, $processed_loh;
  }
  $input->{normal_bam} = $self->normal_bam;
  $input->{tumor_bam} = $self->tumor_bam;
}

sub _run_format_workflow {
  my $self = shift;
  my $input = shift;

  my $lsf_queue = $ENV{'GENOME_LSF_QUEUE_BUILD_WORKER_ALT'};

  my $w = Workflow::Operation->create(
    name => "Format Varscan SNVs",
    operation_type => Workflow::OperationType::Command->get(
      'Genome::Model::Tools::Capture::FormatSnvs')
  );
  $w->parallel_by('variants_file');
  $w->operation_type->lsf_queue($lsf_queue);

  my $log_dir = $self->outdir;
  $w->log_dir($log_dir);

  # Validate the workflow
  my @errors = $w->validate;
  if (@errors) {
    $self->error_message(@errors);
    die "Errors validating workflow\n";
  }

  # Launch workflow
  $self->status_message("Launching workflow now.");
  my $result = Workflow::Simple::run_workflow_lsf(
    $w,
    'variants_file' => $input->{variants_file},
    'outdir' => $self->outdir,
    'append_line' => 1,
  );

  # Collect and analyze results
  unless($result){
    foreach my $error (@Workflow::Simple::ERROR){
      print STDERR $error->error ."\n";
    }
    die $self->error_message("Workflow did not return correctly.");
  }
}

sub _run_process_workflow {
  my $self = shift;
  my $input = shift;

  my $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};

  my $w = Workflow::Operation->create(
    name => "Process Varscan SNVs",
    operation_type => Workflow::OperationType::Command->get(
      'Genome::Model::Tools::Varscan::ProcessSomatic')
  );
  $w->parallel_by('status_file');
  $w->operation_type->lsf_queue($lsf_queue);

  my $log_dir = $self->outdir;
  $w->log_dir($log_dir);

  # Validate the workflow
  my @errors = $w->validate;
  if (@errors) {
    $self->error_message(@errors);
    die "Errors validating workflow\n";
  }

  # Launch workflow
  $self->status_message("Launching workflow now.");
  my $result = Workflow::Simple::run_workflow_lsf(
    $w,
    'status_file' => $input->{format_output},
  );

  #Collect and analyze results
  unless($result){
    foreach my $error (@Workflow::Simple::ERROR){
      print STDERR $error->error ."\n";
    }
    die $self->error_message("Workflow did not return correctly.");
  }
}

sub _run_filter_workflow {
  my $self = shift;
  my $processed_input = shift;
  my $bam = shift;

  my $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};

  my $w = Workflow::Operation->create(
    name => "Filter Varscan SNVs",
    operation_type => Workflow::OperationType::Command->get(
      'Genome::Model::Tools::Somatic::FilterFalsePositives')
  );
  $w->parallel_by('variant_file');
  $w->operation_type->lsf_queue($lsf_queue);

  my $log_dir = $self->outdir;
  $w->log_dir($log_dir);

  # Validate the workflow
  my @errors = $w->validate;
  if (@errors) {
    $self->error_message(@errors);
    die "Errors validating workflow\n";
  }

  # Launch workflow
  $self->status_message("Launching filter workflow now.");
  my $result = Workflow::Simple::run_workflow_lsf(
    $w,
    'variant_file' => $processed_input,
    'bam_file' => $bam,
    'outdir' => $self->outdir,
    'reference' => $self->reference,
    'bam_readcount_version' => "0.6",
  );

  #Collect and analyze results
  unless($result){
    foreach my $error (@Workflow::Simple::ERROR) {
      $self->error_message($error->stdout());
      $self->error_message($error->stderr());
    }
    die $self->error_message("Workflow did not return correctly.");
  }
}

sub execute {
  my $self = shift;
  my %input;
  $self->map_workflow_inputs(\%input);
  $self->_run_format_workflow(\%input);
  $self->_run_process_workflow(\%input);
  $self->_run_filter_workflow($input{processed_somatic}, $input{tumor_bam});
  $self->_run_filter_workflow($input{processed_germline}, $input{tumor_bam});
  $self->_run_filter_workflow($input{processed_loh}, $input{normal_bam});
  return 1;
}

1;

