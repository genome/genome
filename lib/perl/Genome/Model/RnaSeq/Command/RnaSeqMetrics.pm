package Genome::Model::RnaSeq::Command::RnaSeqMetrics;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::RnaSeq::Command::RnaSeqMetrics {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'RNAseq models to generate expression matrix.',
        },
        metrics_tsv_file => {
            doc => '',
            default_value => 'RnaSeqMetrics.tsv',
        },
    ],
    has_optional => [
        normalized_transcript_coverage_file => {
            doc => '',
            default_value => 'NormalizedTranscriptCoverage.tsv',
        },
    ],

};

sub help_synopsis {
    return <<"EOS"
    genome model rna-seq rna-seq-metrics --metrics-tsv-file=FILE MODELS
EOS
}

sub help_brief {
    return "Accumulate RNAseq metrics for models.";
}

sub help_detail {
    return <<EOS
SOMETHING ELSE.
EOS
}

sub execute {
    my $self = shift;
    my @models = $self->models;
    my @non_rna_models = grep { !$_->isa('Genome::Model::RnaSeq') } @models;
    if (@non_rna_models) {
        die('Found a non-RNAseq model: '. Data::Dumper::Dumper(@non_rna_models));
    }
    my @builds;
    my $annotation_build;
    my $reference_build;
    my %model_metrics;
    my @metric_headers = qw/LABEL TOTAL_READS TOTAL_READS_MAPPED TOTAL_READS_UNMAPPED PCT_READS_MAPPED/;
    my @model_metric_keys;
    for my $model (@models) {
        if ( defined($model_metrics{$model->name}) ) {
            die('Multiple models with name: '. $model->name);
        }
        my @model_builds = reverse($model->sorted_builds);
        my $build;
        my $metrics_directory;
        my $metrics_file;
        for (my $i = 0; $i < scalar(@model_builds); $i++) {
            $build = $model_builds[$i];
            $metrics_directory = $build->data_directory .'/metrics';
            $metrics_file = $metrics_directory .'/PicardRnaSeqMetrics.txt';
            if (-s $metrics_file) {
                last;
            }
        }
        push @builds, $build;
        my $model_reference_sequence_build = $model->reference_sequence_build;
        if ($reference_build) {
            unless ($reference_build->id eq $model_reference_sequence_build->id) {
                $self->error_message('Mis-match reference sequence builds!');
                die($self->error_message);
            }
        } else {
            $reference_build = $model_reference_sequence_build;
        }
        my $model_annotation_build = $model->annotation_build;
        if ($annotation_build) {
            unless ($annotation_build->id eq $model_annotation_build->id) {
                $self->error_message('Mis-match annotation builds!');
                die($self->error_message);
            }
        } else {
            $annotation_build = $model_annotation_build;
        }
        unless (-d $metrics_directory) {
            $self->error_message('Missing metrics directory: '. $metrics_directory);
            die($self->error_message);
        }
        unless (-e $metrics_file) {
            $self->error_message('Missing Picard RNAseq metrics file: '. $metrics_file);
            die($self->error_message);
        }
        my $metrics = Genome::Model::Tools::Picard::CollectRnaSeqMetrics->parse_file_into_metrics_hashref($metrics_file);
        unless ($metrics) {
            die('Failed to parse metrics file: '. $metrics_file);
        }
        if ( !@model_metric_keys ) {
            @model_metric_keys = sort keys %{$metrics};
            push @metric_headers, @model_metric_keys;
        } else {
            #TODO: Check that all metrics files have the same headers...
        }
        $model_metrics{$model->name}{metrics} = $metrics;
        my $histo = Genome::Model::Tools::Picard::CollectRnaSeqMetrics->parse_metrics_file_into_histogram_hashref($metrics_file);
        # This is only available in picard v1.52 or greater
        if ($histo) {
            $model_metrics{$model->name}{histogram} = $histo;
        }
    }

    my $metrics_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->metrics_tsv_file,
        separator => "\t",
        headers => \@metric_headers,
    );

    my %transcript_coverage;
    for my $build (@builds) {
        my $name = $build->model->name;
        my $metrics = $model_metrics{$name}{'metrics'};

        my $tophat_stats =  $build->alignment_stats_file;
        unless (-s $tophat_stats) {
            warn('Missing TopHat alignment stats file: '. $tophat_stats);
            next;
        }
        my $tophat_metrics = Genome::Model::Tools::BioSamtools::Tophat2AlignmentStats->parse_alignment_stats_summary_hash_ref($tophat_stats);
        unless (defined($tophat_metrics) || defined($tophat_metrics->{'Total_Reads'})) {
            die('Metrics not parsed correctly from file: '. $tophat_stats);
        }
        my $total_reads = $tophat_metrics->{'TOTAL_READS'};
        if ($total_reads == 0) {
            warn('No reads in stats file: '. $tophat_stats .' or parsed incorrectly : '. Data::Dumper::Dumper($tophat_metrics));
        }
        my $total_reads_mapped = $tophat_metrics->{'TOTAL_READS_MAPPED'};
        my $pct_reads_mapped;
        if ($total_reads) {
            $pct_reads_mapped = ($total_reads_mapped / $total_reads);
        } else {
            $pct_reads_mapped = 0;
        }
        my %summary = (
            LABEL => $name,
            TOTAL_READS => $total_reads,
            TOTAL_READS_MAPPED => $total_reads_mapped,
            TOTAL_READS_UNMAPPED => ($total_reads - $total_reads_mapped),
            PCT_READS_MAPPED => $pct_reads_mapped,
        );
        for my $header (@metric_headers) {
            # This assumes all other values have been set
            if ( !defined($summary{$header}) ) {
                $summary{$header} = $metrics->{$header};
            }
        }
        $metrics_writer->write_one(\%summary);
        if (defined($model_metrics{$name}{'histogram'}) ) {
            my $histo = $model_metrics{$name}{'histogram'};
            for my $position (keys %{$histo}) {
                $transcript_coverage{$position}{$name} = $histo->{$position}{normalized_coverage};
            }
        }
    }

    if ($self->normalized_transcript_coverage_file) {
        my @models = sort keys %model_metrics;
        my @coverage_headers = ('POSITION',@models);
        my $coverage_writer = Genome::Utility::IO::SeparatedValueWriter->create(
            output => $self->normalized_transcript_coverage_file,
            separator => "\t",
            headers => \@coverage_headers,
        );
        for my $position (sort _by_position keys %transcript_coverage) {
            $position =~ /normalized_position-(\d+)/;
            my $position_number = $1;
            my %data = (
                POSITION => $position_number,
            );
            for my $label (@models) {
                $data{$label} = $transcript_coverage{$position}{$label};
            }
        $coverage_writer->write_one(\%data);
        }
    }
    return 1;
}

sub _by_position {
    $a =~ /normalized_position-(\d+)/;
    my $a_digit = $1;
    $b =~ /normalized_position-(\d+)/;
    my $b_digit = $1;
    return $a_digit <=> $b_digit;
}
