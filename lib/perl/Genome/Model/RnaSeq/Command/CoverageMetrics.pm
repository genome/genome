package Genome::Model::RnaSeq::Command::CoverageMetrics;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::RnaSeq::Command::CoverageMetrics {
    is => 'Command::V2',
    has => [
        models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Model group of RNAseq models to generate expression matrix.',
        },
        coverage_tsv_file => {
            doc => '',
            default_value => 'RnaSeqMetrics.tsv',
        },
    ],

};

sub help_synopsis {
    return <<"EOS"
    genome model rna coverage-metrics
EOS
}

sub help_brief {
    return "Accumulate RNA-seq coverage metrics for models.";
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
    my %model_coverage;
    my @coverage_headers = qw/LABEL TOTAL_READS TOTAL_READS_MAPPED UNIQUELY_MAPPED MULTI_MAPPED TOTAL_READS_UNMAPPED PCT_READS_MAPPED/;
    my @model_coverage_keys;
    for my $model (@models) {
        if ( defined($model_coverage{$model->name}) ) {
            die('Multiple models with name: '. $model->name);
        }
        my @model_builds = reverse($model->sorted_builds);
        my $build;
        my $coverage_directory;
        my $coverage_file;
        for (my $i = 0; $i < scalar(@model_builds); $i++) {
            $build = $model_builds[$i];
            $coverage_directory = $build->data_directory .'/coverage';
            $coverage_file = $coverage_directory .'/annotation_squashed_by_gene_STATS.txt';
            if (-s $coverage_file) {
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
        unless (-d $coverage_directory) {
            $self->error_message('Missing coverage directory: '. $coverage_directory);
            die($self->error_message);
        }
        unless (-e $coverage_file) {
            $self->error_message('Missing RefCov RNAseq coverage file: '. $coverage_file);
            die($self->error_message);
        }
        my $coverage_reader = Genome::Utility::IO::SeparatedValueReader->new(input => $coverage_file, separator=> "\t",);
        unless ($coverage_reader) {
            die('Failed to parse coverage file: '. $coverage_file);
        }
        # Presumes these are one line files
        my $coverage = $coverage_reader->next;
        if ( !@model_coverage_keys ) {
            @model_coverage_keys = sort keys %{$coverage};
            push @coverage_headers, @model_coverage_keys;
        } else {
            #TODO: Check that all metrics files have the same headers...
        }
        $model_coverage{$model->name}{coverage} = $coverage;
    }

    my $coverage_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->coverage_tsv_file,
        separator => "\t",
        headers => \@coverage_headers,
    );

    for my $build (@builds) {
        my $name = $build->model->name;
        my $coverage = $model_coverage{$name}{'coverage'};

        my $tophat_stats =  $build->alignment_stats_file;
        my $tophat_fh = Genome::Sys->open_file_for_reading($tophat_stats);
        my %tophat_metrics;
        while (my $line = $tophat_fh->getline) {
            if ($line =~ /^##(.+):\s+(\d+)$/) {
                my $key = $1;
                my $value = $2;
                $key =~ s/ /_/g;
                $tophat_metrics{uc($key)} = $value;
            }
        }
        unless (defined($tophat_metrics{'TOTAL_READS'})) {
            die('Metrics not parsed correctly: '. Data::Dumper::Dumper(%tophat_metrics));
        }
        my %summary = (
            LABEL => $name,
            TOTAL_READS => $tophat_metrics{'TOTAL_READS'},
            TOTAL_READS_MAPPED => $tophat_metrics{'TOTAL_READS_MAPPED'},
            UNIQUELY_MAPPED => $tophat_metrics{'UNIQUE_ALIGNMENTS'},
            MULTI_MAPPED => $tophat_metrics{'MULTIPLE_HIT_READS'},
            TOTAL_READS_UNMAPPED => ($tophat_metrics{'TOTAL_READS'} - $tophat_metrics{'TOTAL_READS_MAPPED'}),
            PCT_READS_MAPPED => ($tophat_metrics{'TOTAL_READS_MAPPED'} / $tophat_metrics{'TOTAL_READS'}),
        );
        for my $header (@coverage_headers) {
            # This assumes all other values have been set
            if ( !defined($summary{$header}) ) {
                $summary{$header} = $coverage->{$header};
            }
        }
        $coverage_writer->write_one(\%summary);
    }

    return 1;
}
