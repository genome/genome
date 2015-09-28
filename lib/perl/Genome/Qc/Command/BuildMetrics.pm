package Genome::Qc::Command::BuildMetrics;

use strict;
use warnings;

use Genome;
use YAML::Syck;
use Test::More qw/no_plan/;

class Genome::Qc::Command::BuildMetrics {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The builds to report QC metrics for.',
        },
        output_file => {
            is => 'Text',
            doc => 'The file path to output build QC metrics.',
        },
    ],
};
sub help_brief {
    'A command to print the QC result metrics for input builds.'
}

sub help_synopsis {
    'Dump QC result metrics from the database to tab-delimited output.'
}

sub help_detail{
    return <<"EOS"
The QC framework stores result metrics in the database for each QC result.  This tool will dump the instrument data QC result metrics for all input builds.  A tab-delimited output of all QC metrics along with build id and instrument data id are output to the terminal.
EOS
}

sub execute {
    my $self = shift;

    my @metrics;

    for my $build ($self->builds) {
        push @metrics, $self->metrics_for_build($build);
    }

    unless (@metrics) {
        $self->error_message('Failed to find QC results for builds!');
        die($self->error_message);
    }

    DumpFile($self->output_file,@metrics);

    return 1;
}

sub metrics_for_build {
    my $self = shift;
    my $build = shift;

    my @metrics;
    my @build_instrument_data = $build->instrument_data;
    my @build_instrument_data_ids = sort(map {$_->id} @build_instrument_data);
    my @qc_results = grep {$_->isa('Genome::Qc::Result')} $build->results;
    for my $qc_result (@qc_results) {
        my $as = $qc_result->alignment_result;
        my @result_instrument_data = $as->instrument_data;
        my @result_instrument_data_ids = sort(map {$_->id} @result_instrument_data);
        # There has to be a way to run diagnostics off
        unless (is_deeply(\@build_instrument_data_ids,\@result_instrument_data_ids)) {
            $self->error_message('Build instrument data and QC result instrument data are not the same!');
            die($self->error_message);
        }
        my %result_metrics = $qc_result->get_unflattened_metrics;
        $result_metrics{instrument_data_count} = scalar(@result_instrument_data_ids);
        $result_metrics{instrument_data_ids} = join(',',@result_instrument_data_ids);
        push @metrics, \%result_metrics;
    }

    return @metrics;
}


1;
