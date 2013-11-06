package Genome::Model::RnaSeq::Command::SpliceJunctionMetrics;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::RnaSeq::Command::SpliceJunctionMetrics {
    is => 'Genome::Command::Base',
    has => [
        models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'Model group of RNAseq models to generate expression matrix.',
        },
        splice_junction_metrics_tsv_file => {
            doc => '',
            default_value => 'SpliceJunctionMetrics.tsv',
        },
    ],

};

sub help_synopsis {
    return <<"EOS"
    genome model rna splice-junction-metrics
EOS
}

sub help_brief {
    return "Accumulate RNA-seq splice junction metrics for models.";
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
    my @junction_headers = qw/LABEL TOTAL_READS TOTAL_READS_MAPPED UNIQUELY_MAPPED MULTI_MAPPED TOTAL_READS_UNMAPPED PCT_READS_MAPPED/;
    my @model_junction_keys;
    for my $model (@models) {
        if ( defined($model_metrics{$model->name}) ) {
            die('Multiple models with name: '. $model->name);
        }
        my @model_builds = reverse($model->sorted_builds);
        my $build;
        my $junctions_directory;
        my $junctions_file;
        for (my $i = 0; $i < scalar(@model_builds); $i++) {
            $build = $model_builds[$i];
            $junctions_directory = $build->junctions_directory;
            unless (-d $junctions_directory) {
                $build = undef;
                next;
            }
            $junctions_file = $junctions_directory .'/summary/Stats.tsv';
            unless (-s $junctions_file) {
                $build = undef;
                next;
            } else {
                last;
            }
        }
        unless ($build) {
            $self->warning_message('Failed to find a build for model \''. $model->id .'\' with splice junction output directory. SKIPPING!');
            next;
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
        unless (-d $junctions_directory) {
            $self->error_message('Missing junctions directory: '. $junctions_directory);
            die($self->error_message);
        }
        unless (-e $junctions_file) {
            $self->error_message('Missing RefCov RNAseq splice junction summary file: '. $junctions_file);
            die($self->error_message);
        }
        my $junction_metrics = Genome::Model::Tools::Transcriptome::SpliceJunctionSummary->parse_summary_file_into_metrics_hash_ref($junctions_file);
        unless (@model_junction_keys) {
            @model_junction_keys = sort keys %{$junction_metrics};
            push @junction_headers, @model_junction_keys;
        }
        $model_metrics{$model->name}{splice_junctions} = $junction_metrics;
    }

    my $splice_junction_metrics_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $self->splice_junction_metrics_tsv_file,
        separator => "\t",
        headers => \@junction_headers,
    );

    for my $build (@builds) {
        my $name = $build->model->name;
        my $splice_junctions = $model_metrics{$name}{'splice_junctions'};

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
        for my $header (@junction_headers) {
            # This assumes all other values have been set
            if ( !defined($summary{$header}) ) {
                $summary{$header} = $splice_junctions->{$header};
            }
        }
        $splice_junction_metrics_writer->write_one(\%summary);
    }

    return 1;
}
