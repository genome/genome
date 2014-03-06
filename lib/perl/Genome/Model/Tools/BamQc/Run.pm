package Genome::Model::Tools::BamQc::Run;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;
use File::Basename qw(basename);

my $DEFAULT_PICARD_VERSION    = Genome::Model::Tools::Picard->default_picard_version;
my $DEFAULT_SAMSTAT_VERSION   = Genome::Model::Tools::SamStat::Base->default_samstat_version;
my $DEFAULT_FASTQC_VERSION    = Genome::Model::Tools::Fastqc->default_fastqc_version; 
my $DEFAULT_ERRORRATE_VERSION = Genome::Model::Tools::BioSamtools::ErrorRate->default_errorrate_version;
my $DEFAULT_SAMTOOLS_VERSION  = Genome::Model::Tools::Sam->default_samtools_version;
my $DEFAULT_SAMTOOLS_MAXIMUM_MEMORY = Genome::Model::Tools::Sam->default_samtools_maximum_memory;

class Genome::Model::Tools::BamQc::Run {
    is => ['Genome::Model::Tools::BamQc::Base'],
    has_input => [
        bam_file => {
            is  => 'Text',
            doc => 'The input BAM file to generate metrics for.'
        },
        output_directory => {
            is  => 'Text',
            doc => 'The directory to write output summary metrics files to.',
        },
        log_directory => {
            is  => 'Text',
            doc => 'The directory to write output summary metrics files to.',
            is_optional => 1,
        },
        reference_sequence => {
            is => 'Text',
            is_optional => 1,
            doc => 'The reference sequence for which the BAM file was aligned to.',
            # GRCh37-lite-build37
            example_values => ['/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fasta',],
        },
        # TODO: Add option for aligner, bwasw would require some additional tools to determine unique alignments
    ],
    has_optional_input => [
        roi_file_path => {
            is  => 'Text',
            doc => 'If supplied ref-cov will run on the supplied regions of interest.',
        },
        roi_file_format => {
            is  => 'Text',
            doc => 'The file format of the supplied ROI',
            valid_values => ['bam','bed'],
        },
        picard_version => {
            is  => 'Text',
            doc => 'The version of Picard to use.  See gmt picard --help for details.',
            default_value => $DEFAULT_PICARD_VERSION,
        },
        samstat_version => {
            is  => 'Text',
            doc => 'The version of SamStat to use. See gmt sam-stat --help for details.',
            default_value => $DEFAULT_SAMSTAT_VERSION,
        },
        fastqc_version => {
            is  => 'Text',
            doc => 'The version of FastQC to use.  See gmt fastqc --help for details.',
            default_value => $DEFAULT_FASTQC_VERSION,
        },
        samtools_version => {
            is  => 'Text',
            doc => 'The version of Samtools to use. See gmt sam --help for details.',
            default_value => $DEFAULT_SAMTOOLS_VERSION,
        },
        samtools_maximum_memory => {
            is  => 'Text',
            doc => 'The maximum amount of memory for Samtools to use. See gmt sam --help for details.',
            default_value => $DEFAULT_SAMTOOLS_MAXIMUM_MEMORY,
        },
        picard_maximum_memory => {
            is  => 'Text',
            doc => 'The maximum amount of memory for Picard to use.  See gmt picard --help for details.',
            default_value => '14',
        },
        picard_maximum_permgen_memory => {
            is  => 'Text',
            doc => 'The maximum amount of permgen memory for Picard to use.  See gmt picard --help for details.',
            default_value => '256',
        },
        picard_max_records_in_ram => {
            is  => 'Text',
            doc => 'The maximum number of records for Picard to keep in memory.  See gmt picard --help for details.',
            default_value => '5000000',
        },
        error_rate => {
            is  => 'Boolean',
            doc => 'A flag to run error rate calculations.',
            default_value => 1,
        },
        error_rate_version => {
            is  => 'Text',
            doc => 'The version of error rate C tool to use.',
            default_value => $DEFAULT_ERRORRATE_VERSION,
        },
        read_length => {
            is  => 'Boolean',
            doc => 'A flag to run read length distribution. This step could take long time on big bam.',
            default_value => 1,
        },
        bam_link_name => {
            is  => 'Text',
            doc => 'The bam symlink name.',
        },
        samstat  => {
            is  => 'Boolean',
            doc => 'A flag to run samstat html-report. Skip this for bwamem/sw bam',
            default_value => 1,
        },
    ],
    has_optional_output => [
        output_metrics_hash_ref => {},
    ],
};

sub execute {
    my $self = shift;
    my %data;

    my $output_directory = $self->output_directory;
    unless (-d $output_directory) {
        unless (Genome::Sys->create_directory($output_directory)) {
            die('Failed to create output directory: '. $output_directory);
        }
    }
    my $log_directory = $self->log_directory || $output_directory;
    unless (-d $log_directory) {
        unless (Genome::Sys->create_directory($log_directory)) {
            die('Failed to create output directory: '. $log_directory);
        }
    }
    unless (-e $self->bam_file) {
        die('Failed to find BAM file: '. $self->bam_file);
    }

    #TODO: Validate BAM sort order or at least throw a warning if not found
    #BAM header tag under @HD SO, valid values: unknown (default), unsorted, queryname and coordinate
    
    my ($bam_basename, $bam_dirname, $bam_suffix) = File::Basename::fileparse($self->bam_file, qw/\.bam/);
    $bam_basename = $self->bam_link_name if $self->bam_link_name;
    my $file_basename = $output_directory .'/'. $bam_basename;

    my $bam_path = $file_basename .'.bam';
    unless (-l $bam_path || -e $bam_path) {
        unless (symlink($self->bam_file,$bam_path)) {
            die('Failed to create symlink '. $bam_path .' -> '. $self->bam_file);
        }
    }

    # SAMTOOLS
    my $bai_path = $file_basename .'.bam.bai';
    unless (-l $bai_path || -e $bai_path) {
        my $bai_file = $self->bam_file .'.bai';
        if (-e $bai_file) {
            unless (symlink($bai_file,$bai_path)) {
                die('Failed to create symlink '. $bai_path .' -> '. $bai_file);
            }
        } 
        else {
            # TODO: test if BAM file is sorted before indexing
            unless (Genome::Model::Tools::Picard::BuildBamIndex->execute(
                use_version            => $self->picard_version,
                maximum_permgen_memory => $self->picard_maximum_permgen_memory,
                maximum_memory         => $self->picard_maximum_memory,
                input_file             => $bam_path,
                output_file            => $bai_path,
            )) {
                die('Failed to index BAM file: '. $bam_path);
            }
        }
    }

    # PICARD MARK DUPLICATES
    my @mrkdup_files = glob($bam_dirname .'/*.metrics');
    unless (@mrkdup_files) {
        # TODO: run MarkDuplicates passing the mrkdup bam file as input to the downstream steps in workflow
        # Not sure because some BAMs are intentionally left unmarked (ex. RNASeq)....
    } 
    else {
        for my $mrkdup_file (@mrkdup_files) {
            my ($mrkdup_basename, $mrkdup_dirname, $mrkdup_suffix) = File::Basename::fileparse($mrkdup_file,qw/\.metrics/);
            my $mrkdup_path = $output_directory .'/'. $mrkdup_basename .'.metrics';
            unless (symlink($mrkdup_file,$mrkdup_path)) {
                die('Failed to create symlink '. $mrkdup_path .' -> '. $mrkdup_file);
            }
        }
    }

    # PICARD METRICS
    my $picard_metrics_basename = $file_basename .'-PicardMetrics';

    # PICARD GC
    my $picard_gc_metrics_file  = $file_basename .'-PicardGC_metrics.txt';
    my $picard_gc_chart_file    = $file_basename .'-PicardGC_chart.pdf';
    my $picard_gc_summary_file  = $file_basename .'-PicardGC_summary.txt';
    
    #ASSUME_SORTED only applied to picard 1.77 and later
    my $picard_gc_assume_sorted = $self->picard_version >= 1.77 ? 1 : 0;

    my %workflow_params = (
        picard_version                 => $self->picard_version,
        output_directory               => $output_directory,
        bam_file                       => $bam_path,
        clean_bam                      => 'none',
        picard_metrics_output_basename => $picard_metrics_basename,
        picard_maximum_permgen_memory  => $self->picard_maximum_permgen_memory,
        picard_maximum_memory          => $self->picard_maximum_memory,
        picard_max_records_in_ram      => $self->picard_max_records_in_ram,
        picard_gc_metrics_file         => $picard_gc_metrics_file,
        picard_gc_chart_file           => $picard_gc_chart_file,
        picard_gc_summary_file         => $picard_gc_summary_file,
        fastqc_version                 => $self->fastqc_version,
    );

    $workflow_params{picard_gc_assume_sorted} = 1 if $picard_gc_assume_sorted; 
    my @output_properties = qw(picard_metrics_result fastqc_result);

    my $flagstat_path = $file_basename .'.bam.flagstat';
    unless (-e $flagstat_path) {
        my $flagstat_file = $self->bam_file .'.flagstat';
        if (-e $flagstat_file) {
            unless (symlink($flagstat_file,$flagstat_path)) {
                die('Failed to create symlink '. $flagstat_path .' -> '. $flagstat_file);
            }
        } 
        else {
            $workflow_params{samtools_version}        = $self->samtools_version;
            $workflow_params{samtools_maximum_memory} = $self->samtools_maximum_memory;
            $workflow_params{samtools_flagstat_path}  = $flagstat_path;
            push @output_properties, 'flagstat_result';
        }
    }

    if ($self->reference_sequence) {
        $workflow_params{reference_sequence} = $self->reference_sequence;
        push @output_properties, 'picard_gc_bias_result';
    }

    if ($self->samstat) {
        $workflow_params{samstat_version} = $self->samstat_version;
        push @output_properties, 'samstat_result';
    }

    if ($self->error_rate) {
        my $error_rate_file = $file_basename .'-ErrorRate.tsv';
        if (-e $error_rate_file) {
            die('Error rate file already exists at: '. $error_rate_file);
        }
        $workflow_params{error_rate_file}    = $error_rate_file;
        $workflow_params{error_rate_version} = $self->error_rate_version;
        push @output_properties, 'error_rate_result';
    }

    if ($self->read_length) {
        $workflow_params{read_length_summary}        = $file_basename .'-ReadLengthSummary.tsv';
        $workflow_params{read_length_histogram}      = $file_basename .'-ReadLengthDistribution.tsv';
        $workflow_params{alignment_length_histogram} = $file_basename .'-AlignmentLengthDistribution.tsv';
        push @output_properties, 'read_length_result';
    }

    if ($self->roi_file_path && $self->reference_sequence) {
        my $refcov_stats_file = $file_basename .'-RefCov_STATS.tsv';
        $workflow_params{roi_file_path}   = $self->roi_file_path;
        $workflow_params{roi_file_format} = $self->roi_file_format;
        $workflow_params{refcov_stats_file} = $refcov_stats_file;
        $workflow_params{refcov_print_headers} = 1;
        $workflow_params{refcov_print_min_max} = 1;

        push @output_properties, 'refcov_result';
    }

    my @input_properties = keys %workflow_params;

    my $workflow = Workflow::Model->create(
        name => 'BamQc',
        input_properties  => \@input_properties,
        output_properties => \@output_properties,
    );

    $workflow->log_dir($log_directory);

    # Samtools flagstat
    if ($workflow_params{samtools_flagstat_path}) {
        my %samtools_flagstat_operation_params = (
            workflow   => $workflow,
            name       => 'Samtools Flagstat',
            class_name => 'Genome::Model::Tools::Sam::Flagstat',
            input_properties => {
                'bam_file'                => 'bam_file',
                'samtools_flagstat_path'  => 'output_file',
                'samtools_version'        => 'use_version',
                'samtools_maximum_memory' => 'maximum_memory',
            },
            output_properties => {
                'result' => 'flagstat_result',
            },
        );
        my $samtools_flagstat_operation = $self->setup_workflow_operation(%samtools_flagstat_operation_params);
    }
    
    # PicardMetrics
    my %picard_metrics_operation_params = (
        workflow   => $workflow,
        name       => 'Collect Picard Metrics',
        class_name => 'Genome::Model::Tools::Picard::CollectMultipleMetrics',
        input_properties => {
            'bam_file'                       => 'input_file',
            'picard_metrics_output_basename' => 'output_basename',
            'picard_version'                 => 'use_version',
            'picard_maximum_memory'          => 'maximum_memory',
            'picard_maximum_permgen_memory'  => 'maximum_permgen_memory',
        },
        output_properties => {
            'result' => 'picard_metrics_result',
        },
    );
    if ($self->reference_sequence) {
        $picard_metrics_operation_params{input_properties}{'reference_sequence'} = 'reference_sequence';
    }
    my $picard_metrics_operation = $self->setup_workflow_operation(%picard_metrics_operation_params);
    my $max_memory = $self->picard_maximum_memory + 2;
    $picard_metrics_operation->operation_type->lsf_resource('-M '. $max_memory .'000000 -R \'select[type==LINUX64 && model!=Opteron250 && tmp>1000 && mem>'. $max_memory.'000] rusage[tmp=1000, mem='. $max_memory.'000]\'');
    
    # PicardGcBias
    if ($workflow_params{reference_sequence}) {
        my %picard_gc_bias_operation_params = (
            workflow   => $workflow,
            name       => 'Collect Picard G+C Bias',
            class_name => 'Genome::Model::Tools::Picard::CollectGcBiasMetrics',
            input_properties => {
                bam_file                      => 'input_file',
                reference_sequence            => 'refseq_file',
                clean_bam                     => 'clean_bam',
                picard_version                => 'use_version',
                picard_maximum_memory         => 'maximum_memory',
                picard_maximum_permgen_memory => 'maximum_permgen_memory',
                picard_max_records_in_ram     => 'max_records_in_ram',
                picard_gc_metrics_file        => 'output_file',
                picard_gc_chart_file          => 'chart_output',
                picard_gc_summary_file        => 'summary_output',
            },
            output_properties => {
                result => 'picard_gc_bias_result',
            },
        );
        $picard_gc_bias_operation_params{input_properties}->{picard_gc_assume_sorted} = 'assume_sorted' 
            if $picard_gc_assume_sorted; 
        
        my $picard_gc_bias_operation = $self->setup_workflow_operation(%picard_gc_bias_operation_params);
        $picard_gc_bias_operation->operation_type->lsf_resource('-M '. $max_memory .'000000 -R \'select[type==LINUX64 && model!=Opteron250 && tmp>1000 && mem>'. $max_memory.'000] rusage[tmp=1000, mem='. $max_memory.'000]\'');
    }
    
    # SamStat
    if ($self->samstat) {
        my %samstat_operation_params = (
            workflow   => $workflow,
            name       => 'SamStat Html Report',
            class_name => 'Genome::Model::Tools::SamStat::HtmlReport',
            input_properties => {
                bam_file        => 'input_files',
                samstat_version => 'use_version',
            },
            output_properties => {
                result => 'samstat_result',
            },
        );
        $self->setup_workflow_operation(%samstat_operation_params);
    }

    # FastQC
    my %fastqc_operation_params = (
        workflow   => $workflow,
        name       => 'FastQC Generate Reports',
        class_name => 'Genome::Model::Tools::Fastqc::GenerateReports',
        input_properties => {
            bam_file         => 'input_files',
            fastqc_version   => 'use_version',
            output_directory => 'report_directory',
        },
        output_properties => {
            result => 'fastqc_result',
        },
    );
    $self->setup_workflow_operation(%fastqc_operation_params);

    # ErrorRate
    if ($self->error_rate) {
        my %error_rate_operation_params = (
            workflow   => $workflow,
            name       => 'BioSamtools Error Rate',
            class_name => 'Genome::Model::Tools::BioSamtools::ErrorRate',
            input_properties => {
                bam_file           => 'bam_file',
                error_rate_file    => 'output_file',
                error_rate_version => 'version',
            },
            output_properties => {
                result => 'error_rate_result',
            },
        );
        my $error_rate_operation = $self->setup_workflow_operation(%error_rate_operation_params);
        $error_rate_operation->operation_type->lsf_resource('-M 8000000 -R \'select[type==LINUX64 && model!=Opteron250 && tmp>1000 && mem>8000] rusage[tmp=1000, mem=8000]\'');
    }
    
    # Read Length
    if ($self->read_length) {
        my %read_length_operation_params = (
            workflow   => $workflow,
            name       => 'BioSamtools Read Length Distribution',
            class_name => 'Genome::Model::Tools::BioSamtools::ReadLengthDistribution',
            input_properties => {
                bam_file                   => 'bam_file',
                read_length_summary        => 'summary_file',
                read_length_histogram      => 'read_length_histogram',
                alignment_length_histogram => 'alignment_length_histogram',
            },
            output_properties => {
                result => 'read_length_result',
            },
        );
        $self->setup_workflow_operation(%read_length_operation_params);
    }

    # RefCov
    # WARNING: MUST USE PERL 5.10.1
    # TODO: The reference_sequence is not a required output
    if ($workflow_params{reference_sequence} && $self->roi_file_path) {
        my %refcov_operation_params = (
            workflow   => $workflow,
            name       => 'RefCov',
            class_name => 'Genome::Model::Tools::RefCov::Standard',
            input_properties => {
                bam_file             => 'alignment_file_path',
                roi_file_path        => 'roi_file_path',
                roi_file_format      => 'roi_file_format',
                refcov_stats_file    => 'stats_file',
                refcov_print_headers => 'print_headers',
                refcov_print_min_max => 'print_min_max',
            },
            output_properties => {
                result => 'refcov_result',
            },
        );
        $self->setup_workflow_operation(%refcov_operation_params);
    }

    my @validation_errors = $workflow->validate;
    unless ($workflow->is_valid) {
        die('Errors encountered while validating workflow: '. join("\n", @validation_errors));
    }
    my $output = Workflow::Simple::run_workflow_lsf($workflow,%workflow_params);
    my @execution_errors = @Workflow::Simple::ERROR;
    if (@execution_errors) {
        for (@execution_errors) {
            print STDERR $_->error ."\n";
        }
        return;
    }

    # SUMMARY
    # TODO: Eventually create and xls spreadsheet or consolidated PDF report
    my $flagstat_data        = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    $data{'FlagstatMetrics'} = $flagstat_data;

    my $insert_size_file = $picard_metrics_basename .'.insert_size_metrics';
    if (-e $insert_size_file) {
        my $insert_size_data  = Genome::Model::Tools::Picard::CollectInsertSizeMetrics->parse_file_into_metrics_hashref($insert_size_file);
        my $insert_size_histo = Genome::Model::Tools::Picard::CollectInsertSizeMetrics->parse_metrics_file_into_histogram_hashref($insert_size_file);
        $data{'InsertSizeMetrics'}   = $insert_size_data;
        $data{'InsertSizeHistogram'} = $insert_size_histo;
    }

    my $alignment_summary_file = $picard_metrics_basename .'.alignment_summary_metrics';
    if (-e $alignment_summary_file) {
        my $alignment_summary_data = Genome::Model::Tools::Picard::CollectAlignmentSummaryMetrics->parse_file_into_metrics_hashref($alignment_summary_file);
        $data{'AlignmentSummaryMetrics'} = $alignment_summary_data;
    }

    my $quality_score_file = $picard_metrics_basename .'.quality_distribution_metrics';
    if (-e $quality_score_file) {
        my $quality_score_histo = Genome::Model::Tools::Picard::QualityScoreDistribution->parse_metrics_file_into_histogram_hashref($quality_score_file);
        $data{'QualityScoreHistogram'} = $quality_score_histo;
    }

    my $quality_cycle_file = $picard_metrics_basename .'.quality_by_cycle_metrics';
    if (-e $quality_cycle_file) {
        my $quality_cycle_histo = Genome::Model::Tools::Picard::MeanQualityByCycle->parse_metrics_file_into_histogram_hashref($quality_cycle_file);
        $data{'MeanQualityByCycleHistogram'} = $quality_cycle_histo;
    }

    if (-e $picard_gc_metrics_file) {
        my $gc_data = Genome::Model::Tools::Picard::CollectGcBiasMetrics->parse_file_into_metrics_hashref($picard_gc_metrics_file);
        $data{'GcBiasMetrics'} = $gc_data;
    }

    if (-e $picard_gc_summary_file) {
        my $gc_data = Genome::Model::Tools::Picard::CollectGcBiasMetrics->parse_file_into_metrics_hashref($picard_gc_summary_file);
        $data{'GcBiasSummary'} = $gc_data;
    }

    if ($self->roi_file_path) {
        my $refcov_stats = Genome::Utility::IO::SeparatedValueReader->create(
            input => $workflow_params{refcov_stats_file},
            separator => "\t",
        );
        while (my $refcov_data = $refcov_stats->next) {
            if ($data{'RefCovMetrics'}{$refcov_data->{'name'}}) {
                die('Multiple RefCov entries found.  Probably from multiple min_depth or wingspan filters.');
            } 
            else {
                $data{'RefCovMetrics'}{$refcov_data->{'name'}} = $refcov_data;
            }
        }
    }

    for my $mrkdup_file (@mrkdup_files) {
        my ($mrkdup_basename,$mrkdup_dir,$mrkdup_suffix) = File::Basename::fileparse($mrkdup_file,qw/\.metrics/);
        my $mrkdup_symlink = $output_directory.'/'. $mrkdup_basename . $mrkdup_suffix;
        unless (-e $mrkdup_symlink) {
            symlink($mrkdup_file,$mrkdup_symlink) || die('Failed to symlink '. $mrkdup_symlink .' -> '. $mrkdup_file);
        }
        my $mrkdup_data = Genome::Model::Tools::Picard::MarkDuplicates->parse_file_into_metrics_hashref($mrkdup_symlink);
        for my $library (keys %{$mrkdup_data}) {
            if (defined($data{'MarkDuplicatesMetrics'}{$library})) {
                die('More than one MarkDuplicates file found for library '. $library .' in '. Data::Dumper::Dumper(@mrkdup_files));
            }
            $data{'MarkDuplicatesMetrics'}{$library} = $mrkdup_data->{$library};
        }
    }

    $self->output_metrics_hash_ref(\%data);
    return 1;
}

sub setup_workflow_operation {
    my ($self, %params) = @_;

    my $workflow   = delete($params{'workflow'});
    my $name       = delete($params{'name'});
    my $class_name = delete($params{'class_name'});
    my $input_properties  = delete($params{'input_properties'});
    my $output_properties = delete($params{'output_properties'});

    my $input_connector  = $workflow->get_input_connector;
    my $output_connector = $workflow->get_output_connector;

    my $operation = $workflow->add_operation(
        name => $name,
        operation_type => Workflow::OperationType::Command->create(
            command_class_name => $class_name,
        )
    );

    for my $left_property (keys %{$input_properties}) {
        my $right_property = $input_properties->{$left_property};
        $workflow->add_link(
            left_operation  => $input_connector,
            left_property   => $left_property,
            right_operation => $operation,
            right_property  => $right_property,
        );
    }
    for my $left_property (keys %{$output_properties}) {
        my $right_property = $output_properties->{$left_property};
        $workflow->add_link(
            left_operation  => $operation,
            left_property   => $left_property,
            right_operation => $output_connector,
            right_property  => $right_property,
        );
    }
    return $operation;
}



1;
