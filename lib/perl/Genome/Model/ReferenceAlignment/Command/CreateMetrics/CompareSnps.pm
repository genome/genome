package Genome::Model::ReferenceAlignment::Command::CreateMetrics::CompareSnps;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::ReferenceAlignment::Command::CreateMetrics::CompareSnps {
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
    doc => 'Creates metrics for a build from the compare_snps file',
};

sub help_brief { return 'Parses the compare_snps file and creates metrics on the build' };
sub help_synopsis { help_brief() };
sub help_detail {
    return 'If the given build in lane qc, metrics are generated based on the compare snps file. If the build ' .
    'is regular reference alignment, then the lane qc builds are found and metrics are created for them.';
}

# Fields in the compare snps file that should be made into metrics
sub compare_snps_fields_to_metrics {
    return qw/
        SnpsCalled
        WithGenotype
        OverallConcord
    /;
}

sub execute {
    my $self = shift;

    my @builds;
    if ($self->build->model->is_lane_qc) {
        push @builds, $self->build;
    }
    else {
        push @builds, map { $_->lane_qc_build } $self->build->instrument_data;
    }

    for my $build (@builds) {
        my $output_file = $build->compare_snps_file;
        confess "No compare_snps file found for build " . $build->id . " at $output_file" unless -e $output_file;

        my @output_columns = Genome::Model::Tools::Analysis::LaneQc::CompareSnps->output_columns;
        my $svr = Genome::Utility::IO::SeparatedValueReader->create(
            headers => \@output_columns,
            separator => "\t",
            is_regex => 1,
            ignore_extra_columns => 0,
            input => $output_file
        );
        confess "Could not create separated value reader for $output_file!" unless $svr;

        $svr->next; # Skip headers
        my $line = $svr->next; # Report only has one line of output that we care about

        for my $field ($self->compare_snps_fields_to_metrics) {
            my $metric_value = $line->{$field};
            next unless defined $metric_value;

            my $metric_name = lc( 'compare_snps_'.Genome::Utility::Text::camel_case_to_string($field, '_') );
            my $metric = Genome::Model::Metric->get(
                build_id => $build->id,
                name => $metric_name,
            );
            if ($metric) {
                $metric->value($metric_value);
            }
            else {
                Genome::Model::Metric->create(
                    build_id => $build->id,
                    name => $metric_name,
                    value => $metric_value,
                );
            }
        }
    }

    return 1;
}

1;

