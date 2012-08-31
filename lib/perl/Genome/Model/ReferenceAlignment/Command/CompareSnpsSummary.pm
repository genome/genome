package Genome::Model::ReferenceAlignment::Command::CompareSnpsSummary;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::ReferenceAlignment::Command::CompareSnpsSummary {
    is => 'Genome::Command::Base',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            id_by => 'build_id',
            shell_args_position => 1,
        },
        build_id => {
            is => 'Number',
        },
    ],
    has_optional => [
        force_new_metrics => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, the compare_snps files are re-parsed and new metrics are made instead of using existing metrics',
        },
    ],
};

sub help_brief { return 'Prints a summary of all the compare_snps metrics pulled from underlying lane qc builds'; }
sub help_synopsis { return help_brief(); }
sub help_detail {
    return 'Looks for compare_snps metrics on lane qc builds associated with the supplied build and generates them if none ' .
    'are found. All the metrics are then pooled together and a report printed';
}

sub summary_headers {
    return qw/
        instrument_data_id
        flow_cell
        qc_build
        lane
        snps_called 
        with_genotype
        overall_concord
    /;
}

sub qc_metrics {
    return qw/
        compare_snps_snps_called
        compare_snps_with_genotype
        compare_snps_overall_concord
    /;
}

sub fields_to_average {
    return qw/
        snps_called
        with_genotype
        overall_concord
    /;
}

sub execute {
    my $self = shift;
    my @instrument_data = $self->build->instrument_data;
    confess 'Found no instrument data assigned to build!' unless @instrument_data;

    my %report_info;
    for my $instrument_data (@instrument_data) {
        $report_info{$instrument_data->id}{'flow_cell'} = $self->get_flow_cell_for_instrument_data($instrument_data);
        $report_info{$instrument_data->id}{'lane'} = $self->get_lane_for_instrument_data($instrument_data);

        my $qc_build = $instrument_data->lane_qc_build;
        if ($qc_build) {
            $report_info{$instrument_data->id}{'qc_build'} = $qc_build->id;
        }
        else {
            $report_info{$instrument_data->id}{'qc_build'} = 'not_found';
            next;
        }

        my @metrics = $self->get_compare_snps_metrics($qc_build);
        next unless @metrics;
        for my $metric (@metrics) {
            my $name = $metric->name;
            $name =~ s/compare_snps_//; # compare_snps is added to metric name to help idenify it, remove it here
            my $value = $metric->value;
            $value =~ s/\%//; # percent sign messes up averaging, remove it
            $report_info{$instrument_data->id}{$name} = $value;
        }
    }

    $self->print_report(%report_info);
    return 1;
}

sub print_report {
    my ($self, %report) = @_;
    my %totals; # Used for averaging later

    print join("\t", $self->summary_headers) . "\n"; 

    for my $instrument_data_id (sort keys %report) {
        my @info;
        push @info, $instrument_data_id;
        for my $field ($self->summary_headers) {
            next if $field eq 'instrument_data_id';
            my $value = $report{$instrument_data_id}{$field};
            $value ||= '-';
            push @info, $value;

            if (grep { $field eq $_ } $self->fields_to_average and $value ne '-') {
                $totals{$field}{'total'} += $value;
                $totals{$field}{'num_values'}++;
            }
        }
        print join("\t", @info) . "\n";
    }
    
    print "\nAverages:\n";
    for my $field (sort keys %totals) {
        my $value = sprintf("%.3f", $totals{$field}{'total'} / $totals{$field}{'num_values'});
        print "$field: $value\n";
    }
    return 1;
}

sub get_lane_for_instrument_data {
    my ($self, $instrument_data) = @_;
    my $lane;
    if ($instrument_data->can('lane') and defined $instrument_data->lane) {
        $lane = $instrument_data->lane;
    }
    return $lane;
}

sub get_flow_cell_for_instrument_data {
    my ($self, $instrument_data) = @_;
    my $flow_cell_id;
    if ($instrument_data->can('flow_cell_id') and defined $instrument_data->flow_cell_id) {
        $flow_cell_id = $instrument_data->flow_cell_id;
    }
    return $flow_cell_id;
}

sub get_compare_snps_metrics {
    my ($self, $qc_build) = @_;
    my @metrics = Genome::Model::Metric->get(
        build_id => $qc_build->id,
        name => [$self->qc_metrics],
    );
    return @metrics if @metrics and not $self->force_new_metrics;
    return $self->create_and_retrieve_qc_metrics($qc_build);
}

sub create_and_retrieve_qc_metrics {
    my ($self, $qc_build) = @_;

    my $rv = Genome::Model::ReferenceAlignment::Command::CreateMetrics::CompareSnps->execute(
        build_id => $qc_build->id,
    );
    return unless $rv;

    return Genome::Model::Metric->get(
        build_id => $qc_build->id,
        name => [$self->qc_metrics],
    );
}

1;

