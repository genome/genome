package Genome::Qc::Command::BuildMetrics;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Qc::Command::BuildMetrics {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The builds to report QC metrics for.',
        },
    ]
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
        my @qc_results = grep {$_->isa('Genome::Qc::Result')} $build->results;
        for my $qc_result (@qc_results) {
            my $as = $qc_result->alignment_result;
            my @instrument_data = $as->instrument_data;
            if (@instrument_data > 1) {
                $self->error_message('Please add support for merged alignment results and multiple instrument data QC results!');
                die($self->error_message);
            }
            my %metrics = $qc_result->get_metrics;
            $metrics{build_id} = $build->id;
            $metrics{instrument_data_id} = $instrument_data[0]->id;
            push @metrics, \%metrics;
        }
    }
    unless (@metrics) {
        $self->error_message('Failed to find QC results for builds!');
        die($self->error_message);
    }
    my @headers = List::MoreUtils::uniq(sort map { keys %$_ } @metrics);
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
    );
    for (@metrics) {
        $writer->write_one($_);
    }
    $writer->output->close;
    
    return 1;
}

1;
