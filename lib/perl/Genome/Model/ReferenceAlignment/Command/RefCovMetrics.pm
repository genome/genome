package Genome::Model::ReferenceAlignment::Command::RefCovMetrics;

use strict;
use warnings 'FATAL';

use Genome;

require List::MoreUtils;
use Params::Validate ':types';

class Genome::Model::ReferenceAlignment::Command::RefCovMetrics {
    is => 'Command::V2',
    has_input => {
        models => {
            is => 'Genome::Model',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'Models to get the reference coverage metrics.',
        },
        type => {
            is => 'Text',
            valid_values => [qw/ alignment coverage /],
            doc => 'Type of metrics to retrieve from the reference coverage result',
        },
    },
    has_optional_input => {
        separator => {
            is => 'Text',
            default_value => "\t",
            doc => 'Separator to use. Default value is tab (\t).',
        },
        row_ids => {
            is => 'Text',
            is_many => 1,
            valid_values => [qw/
                model_id model_name
                result_id
                sample_name sample_common_name
            /],
            default_value => [qw/ model_name /],
            doc => 'The values to use to as the row identifier(s).',
        },
    },
    has_optional_output => {
        output_path => {
            is => 'Text',
            default_value => '-',
            doc => 'Path to send output.',
        },
    },
};

sub help_brief { 'get ref cov metrics from models' }
sub help_detail {
    return <<HELP;
This command will get reference coverage from a model's succeeded builds. It goes directly to the refcof software result, and does not look at the model/build metrics.
HELP
}

sub execute {
    my $self = shift;

    my ($metrics_method, $extract_metrics_method) = $self->_resolve_methods;
    my $fh = Genome::Sys->open_file_for_writing($self->output_path);

    my (@headers, @data);
    MODEL: for my $model ( $self->models ) {
        BUILD: for my $build ( $model->builds(status => 'Succeeded', -order => 'date_completed') ) {
            if ( not $build->is_current ) {
                next BUILD;
            }

            my @results = $build->results;
            next BUILD if not @results;
            my @coverage_stats = List::MoreUtils::uniq(
                grep { $_->isa('Genome::InstrumentData::AlignmentResult::Merged::CoverageStats') } @results
            );
            next BUILD if not @coverage_stats;

            RESULT: for my $result ( @coverage_stats ) {
                my $metrics = $result->$metrics_method;
                next RESULT if not $metrics;
                push @headers, keys %$metrics;
                my $row = $self->$extract_metrics_method($metrics);
                $self->_add_row_ids($metrics, {model => $model, result => $result});
                push @data, $row;
                next MODEL;
            }
        }
        $self->warning_message('No results for model: '.$model->__display_name__);
    }

    return 1 if not @data;

    my @headers = List::MoreUtils::uniq(sort @headers);
    unshift @headers, $self->row_ids;
    $fh->print( join($self->separator, @headers)."\n" );
    for my $data ( @data ) {
        $fh->print( join($self->separator, map { $data->{$_} } @headers)."\n" );
    }

    $fh->close;

    return 1;
}

sub _resolve_methods {
    my $self = shift;
    if ( $self->type eq 'alignment' ) {
        return (qw/ alignment_summary_hash_ref _extract_alignment_metrics /);
    }
    else { 
        return (qw/ coverage_stats_summary_hash_ref _extract_coverage_metrics /);
    }
}

sub _add_row_ids {
    my ($self, $metrics, $objects)= Params::Validate::validate_pos(@_, {type => SCALAR}, {type => HASHREF});

    my @samples = List::MoreUtils::uniq(
        map { $_->sample } $objects->{result}->alignment_result->instrument_data
    );
    $objects->{sample} = ( @samples == 1 )
    ? $samples[0]
    : $objects->{model}->subject; # this could be an individual

    my @ids;
    for my $row_id ( $self->row_ids ) {
        my ($object, $method) = split(/_/, $row_id, 2);
        $metrics->{$row_id} = $objects->{$object}->$method;
    }

    return @ids;
}

sub _extract_alignment_metrics {
    my ($self, $metrics) = @_;

    my %row;
    for my $wingspan ( keys %$metrics ) {
        for my $name ( keys %{$metrics->{$wingspan}} ) {
            my $key = join('_', 'alignment-wingspan', $wingspan, $name);
            $row{$key} = $metrics->{$wingspan}->{$name};
        }
    }

    return \%row;
}

sub _extract_coverage_metrics {
    my ($self, $metrics) = @_;

    my %row;
    for my $wingspan ( keys %$metrics ) {
        for my $minimum_depth ( keys %{$metrics->{$wingspan}} ) {
            for my $name ( keys %{$metrics->{$wingspan}->{$minimum_depth}} ) {
                my $key = join('_', 'coverage-wingspan', $wingspan, $minimum_depth, $name);
                $row{$key} = $metrics->{$wingspan}->{$minimum_depth}->{$name};
            }
        }
    }

    return \%row;
}

1;

