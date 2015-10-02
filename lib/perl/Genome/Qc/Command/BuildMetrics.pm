package Genome::Qc::Command::BuildMetrics;

use strict;
use warnings;

use Genome;
use YAML::Syck qw(DumpFile);

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
            doc => 'The file path to output build QC metrics as YAML.',
        },
    ],
};
sub help_brief {
    'A command to print the QC result metrics for input builds.'
}

sub help_synopsis {
    'Dump QC result metrics from the database to a YAML output file.'
}

sub help_detail{
    return <<"EOS"
The QC framework stores result metrics in the database for each QC result.  This tool will dump the QC result metrics for all input builds.  A YAML format output of all QC metrics along with build id and instrument data ids are output to the defined output file.
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
    my $build_instdata_set = Set::Scalar->new($build->instrument_data);
    my @qc_results = grep {$_->isa('Genome::Qc::Result')} $build->results;
    for my $qc_result (@qc_results) {
        my $as = $qc_result->alignment_result;
        my $result_instdata_set = Set::Scalar->new($as->instrument_data);
        if ($build_instdata_set->is_equal($result_instdata_set)) {
            my %result_metrics = $qc_result->get_unflattened_metrics;
            $result_metrics{build_id} = $build->id;
            $result_metrics{instrument_data_count} = $result_instdata_set->size;
            $result_metrics{instrument_data_ids} = join(',',map {$_->id} $result_instdata_set->members);
            push @metrics, \%result_metrics;
        } else {
            $self->error_message('Build and QC result instrument data are not the same!');
            die($self->error_message);
        }
    }

    return @metrics;
}


1;
