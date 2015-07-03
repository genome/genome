package Genome::Model::ReferenceAlignment::Command::RefCovMetrics;

use strict;
use warnings 'FATAL';

use Genome;

require List::MoreUtils;

class Genome::Model::ReferenceAlignment::Command::RefCovMetrics {
    is => 'Command::V2',
    has_input => {
        models => {
            is => 'Genome::Model::ReferenceAlignment',
            is_many => 1,
            shell_args_position => 1,
            require_user_verify => 0,
            doc => 'Models to get the reference coverage metrics.',
        },
        type => {
            is => 'Text',
            default_value => 'coverage',
            valid_values => [qw/ alignment coverage /],
            doc => 'Type of metrics to retrieve from the reference coverage result',
        },
    },
    has_optional_input => {
        separator => {
            is => 'Text',
            default_value => "\t",
            doc => 'Separator to use.',
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

sub execute {
    my $self = shift;

    my ($metrics_method, $extract_metrics_method) = $self->_resolve_methods;
    my $fh = Genome::Sys->open_file_for_writing($self->output_path);

    my (@data, @headers);
    MODEL: for my $model ( $self->models ) {
        BUILD: for my $build ( $model->builds(status => 'Succeeded', -order => 'date_completed') ) {
            my $metrics = $build->$metrics_method;
            next BUILD if not $metrics;
            my $row = $self->$extract_metrics_method($metrics);
            $row->{model_name} = $model->name;
            push @data, $row;
            push @headers, keys %$row;
            next MODEL;
        }
    }

    my @headers = List::MoreUtils::uniq(sort map { keys %$_ } @data);
    my $idx = List::MoreUtils::firstidx { $_ eq 'model_name' } @headers;
    splice(@headers, $idx, 1);
    unshift @headers, 'model_name';
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

sub _extract_alignment_metrics {
    my ($self, $metrics) = @_;

    my %row;
    for my $wingspan ( keys %$metrics ) {
        print "$wingspan\n";
        for my $name ( keys %{$metrics->{$wingspan}} ) {
            my $key = join('_', 'alignment-wingspan', $wingspan, $name);
            $row{$name} = $metrics->{$wingspan}->{$name};
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

