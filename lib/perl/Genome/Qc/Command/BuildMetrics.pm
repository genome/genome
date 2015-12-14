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

    DumpFile($self->output_file,\@metrics);

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
            $result_metrics{instrument_data_ids} = join(',',sort map {$_->id} $result_instdata_set->members);
            if ($result_metrics{PAIR}) {
                # Calculate Duplication Rate
                if ( defined($result_metrics{reads_marked_duplicates}) ) {
                    $result_metrics{DUPLICATION_RATE} = $result_metrics{'reads_marked_duplicates'}
                        / $result_metrics{PAIR}->{PF_READS_ALIGNED};
                } else {
                    $self->error_message('Missing samtools reads_marked_duplicates!');
                    die($self->error_message);
                }
                # Calculate Haploid Coverage
                $self->_calculate_haploid_coverage($build,\%result_metrics);
            } else {
                $self->error_message('Missing CollectAlignmentSummaryMetrics PAIR category.');
                die($self->error_message);
            }
            push @metrics, \%result_metrics;
        } else {
            $self->error_message('Build and QC result instrument data are not the same!');
            die($self->error_message);
        }
    }

    return @metrics;
}

sub _calculate_haploid_coverage {
    my $self = shift;
    my $build = shift;
    my $result_metrics = shift;

    if ( !defined($result_metrics->{GENOME_TERRITORY}) ) {
        if ( !defined($build->reference_sequence_build->get_metric('GENOME_TERRITORY')) ) {
            my $calc_genome_territory_cmd = Genome::Model::ReferenceSequence::Command::CalculateGenomeTerritory->create(
                reference_sequence_build => $build->reference_sequence_build,
            );
            unless ($calc_genome_territory_cmd->execute) {
                $self->error_message('Failed to execute CalcuateGenomeTerritory command!');
                die($self->error_message);
            }
        }
        $result_metrics->{GENOME_TERRITORY} = $build->reference_sequence_build->get_metric('GENOME_TERRITORY');
    }
    $result_metrics->{HAPLOID_COVERAGE} = (
        $result_metrics->{PAIR}->{PF_ALIGNED_BASES} * ( 1 - $result_metrics->{DUPLICATION_RATE} )
    ) / $result_metrics->{GENOME_TERRITORY};

    return 1;
}

1;
