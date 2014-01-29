package Genome::Model::RnaSeq::Command::CufflinksFpkmTrackingStats;

use strict;
use warnings;

use Genome;
use Statistics::Descriptive;

class Genome::Model::RnaSeq::Command::CufflinksFpkmTrackingStats {
    is => 'Genome::Command::Base',
    has_input => [
        models => {
            is => 'Genome::Model::RnaSeq',
            is_many => 1,
            shell_args_position => 1,
            doc => 'RNAseq models to generate expression matrix.',
        },
        gene_fpkm_tracking_stats_tsv_file => {
            doc => 'The output tsv file of gene-level FPKM tracking stats.',
        },
        isoform_fpkm_tracking_stats_tsv_file => {
            doc => 'The output tsv file of isoform-level FPKM tracking stats.',
            is_optional => 1,
        },
        model_identifier => {
            is_optional => 1,
            default_value => 'name',
            valid_values => ['name','subject_name'],
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
    
EOS
}

sub help_brief {
    return "";
}

sub help_detail {
    return <<EOS

EOS
}


sub execute {
    my $self = shift;

    my @models = $self->models;
    my %feature_types = (
        gene => 'tracking_id',
        isoform => 'tracking_id',
    );
    my @builds;
    my $annotation_build;
    my $reference_build;
    my %model_identifiers;
    my $method = $self->model_identifier;
    for my $model (@models) {
        my $build = $model->last_succeeded_build;
        unless ($build) {
            $build = $model->latest_build;
            unless ($build) {
                die('Failed to find build for model: '. $model->id);
            }
        }
        push @builds, $build;
        my $model_reference_sequence_build = $model->reference_sequence_build;
        if ($reference_build) {
            unless ($reference_build->id eq $model_reference_sequence_build->id) {
                die('Mis-match reference sequence builds!');
            }
        } else {
            $reference_build = $model_reference_sequence_build;
        }
        my $model_annotation_build = $model->annotation_build;
        if ($annotation_build) {
            unless ($annotation_build->id eq $model_annotation_build->id) {
                die('Mis-match annotation builds!');
            }
        } else {
            $annotation_build = $model_annotation_build;
        }
    }
    my @fpkm_tracking_headers;
    my %model_stats;
    for my $build (@builds) {
        my $identifier = $build->model->$method;
        if ( defined($model_stats{$identifier}) ) {
            die('Multiple models with '. $method .' : '. $identifier);
        }
        for my $feature_type (keys %feature_types) {
            my $output_file_method = $feature_type .'_fpkm_tracking_stats_tsv_file';
            # Skip a feature type if the output file method is not defined
            unless  ($self->$output_file_method) { next; }

            my $fpkm_tracking = $build->data_directory .'/expression/'. $feature_type .'s.fpkm_tracking';
            unless (-e $fpkm_tracking) {
                die ('Failed to find '. $feature_type .' FPKM file: '. $fpkm_tracking);
            }
            $self->debug_message('Generating FPKM tracking stats for: '. $fpkm_tracking);
            my $fpkm_tracking_stats_cmd = Genome::Model::Tools::Cufflinks::FpkmTrackingStats->create(
                fpkm_tracking_file => $fpkm_tracking,
            );
            unless ($fpkm_tracking_stats_cmd->execute) {
                die('Failed to execute FpkmTrackingStats for file: '. $fpkm_tracking);
            }
            my $stats = $fpkm_tracking_stats_cmd->_stats_hash_ref;
            unless (@fpkm_tracking_headers) {
                @fpkm_tracking_headers = sort keys %{$stats};
            }
            $model_stats{$identifier}{$feature_type} = $stats;
        }
    }
    $self->debug_message('Printing output FPKM tracking stats file...');
    my @output_headers = ('model_id',@fpkm_tracking_headers);
    for my $feature_type (keys %feature_types) {
        my $output_file_method = $feature_type .'_fpkm_tracking_stats_tsv_file';
        unless ($self->$output_file_method) { next; }
        my $output_file = $self->$output_file_method;
        my $tsv_writer = Genome::Utility::IO::SeparatedValueWriter->create(
            output => $output_file,
            separator => "\t",
            headers => \@output_headers,
        );
        unless ($tsv_writer) {
            die('Failed to open '. $feature_type .' FPKM output file: '. $output_file);
        }
        for my $model_id (sort keys %model_stats) {
            my $stats = $model_stats{$model_id}{$feature_type};
            $stats->{model_id} = $model_id;
            $tsv_writer->write_one($stats);
        }
    }
    $self->debug_message('Finished!');
    return 1;
}
