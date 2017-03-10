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
        output_format => {
            is => 'Text',
            doc => 'The output file format.',
            valid_values => ['yaml','tsv'],
        },
        somatic => {
            is => 'Boolean',
            doc => 'build is somatic or not',
            default_value => 0,
        }
    ],
    has_optional => [
        _tsv_output_profile => {
            is => 'Text',
            doc => 'The output TSV profile used to determine metrics reported in file.',
            valid_values => ['wgs','exome'],
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
        push @metrics, $self->calculate_qc_metrics_for_build($build);
    }

    unless (@metrics) {
        $self->error_message('Failed to find QC results for builds!');
        die($self->error_message);
    }

    if ($self->output_format eq 'yaml') {
        DumpFile($self->output_file,@metrics);
    } elsif ($self->output_format eq 'tsv') {
        $self->output_metrics_as_tsv(\@metrics);
    }

    return 1;
}

sub calculate_qc_metrics_for_build {
    my $self = shift;
    my $build = shift;

    my @metrics;
    my $build_instdata_set = Set::Scalar->new($build->instrument_data);
    my @qc_results = grep {$_->isa('Genome::Qc::Result')} $build->results;
    for my $qc_result (@qc_results) {
        my $config_type = $qc_result->qc_config->type;
        if (!defined($self->_tsv_output_profile)) {
            $self->_tsv_output_profile($config_type);
        } elsif ($self->_tsv_output_profile ne $config_type) {
            $self->fatal_message('Unable to summarize QC metrics for builds with multiple qc config types!');
        }
        
        my $as = $qc_result->alignment_result;
        my $result_instdata_set = Set::Scalar->new($as->instrument_data);
        if ($self->somatic or $build_instdata_set->is_equal($result_instdata_set)) {
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
                # Calculate Haploid Coverage for WGS
                if ($qc_result->qc_config->type eq 'wgs') {
                    $self->_calculate_haploid_coverage($build,\%result_metrics);
                }
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


sub output_metrics_as_tsv {
    my $self = shift;
    my $qc_metric_results = shift;
    
    my @headers = qw/
                       sample_name
                       build_id
                       instrument_data_count
                       instrument_data_ids
                       bam_path
                       pf_reads
                       aligned_bases
                       alignment_rate
                       duplication_rate
                       first_of_pair_mismatch_rate
                       second_of_pair_mismatch_rate
                       microarray_chipmix
                       microarray_freemix
                   /;
    my @exome_headers = qw/
                              mean_target_coverage
                              pct_usable_bases_on_target
                              pct_target_bases_20x
                              pct_exc_off_target
                              pct_exc_dupe
                              pct_exc_overlap
                              pct_exc_mapq
                              pct_exc_baseq
                          /;
    my @wgs_headers = qw/
                             chipmix
                             freemix
                             haploid_coverage
                        /;

    # The results themselves have a type, eg. wgs or exome
    # Rather than relying on a user input profile, could use the result type instead
    my @data;
    if ($self->_tsv_output_profile eq 'wgs') {
        #TODO: Add additional WGS specific metrics?
        push @headers, @wgs_headers;
        for my $qc_metric_result (@$qc_metric_results) {
            my $data = $self->_base_hash_ref_for_qc_metric_result($qc_metric_result);
            $data->{'chipmix'} = $qc_metric_result->{ALL}->{CHIPMIX};
            $data->{'freemix'} = $qc_metric_result->{ALL}->{FREEMIX};
            $data->{'haploid_coverage'} = $qc_metric_result->{HAPLOID_COVERAGE};
            $self->_add_microarray_vbid_result_metrics_to_hash_ref($qc_metric_result,$data);
            push @data, $data;
        }
    } elsif ($self->_tsv_output_profile eq 'exome') {
        push @headers, @exome_headers;
        for my $qc_metric_result(@$qc_metric_results) {
            my $data = $self->_base_hash_ref_for_qc_metric_result($qc_metric_result);
            for my $column (@exome_headers) {
                $data->{$column} = $qc_metric_result->{uc $column};
            }
            $self->_add_microarray_vbid_result_metrics_to_hash_ref($qc_metric_result,$data);
            push @data, $data;
        }
    } else {
        die('Invalid TSV output QC profile!');
    }
    $self->_write_metrics_hash_refs_to_tsv_file(\@headers,\@data);
    return 1;
}

sub _base_hash_ref_for_qc_metric_result {
    my $self = shift;
    my $qc_metric_result = shift;

    my $build = Genome::Model::Build->get($qc_metric_result->{build_id});
    my %data = (
        sample_name => $build->model->subject->name,
        build_id => $qc_metric_result->{build_id},
        instrument_data_count => $qc_metric_result->{instrument_data_count},
        instrument_data_ids => $qc_metric_result->{instrument_data_ids},
        bam_path => $build->merged_alignment_result->bam_path,
        pf_reads => $qc_metric_result->{PAIR}->{PF_READS},
        aligned_bases => $qc_metric_result->{PAIR}->{PF_ALIGNED_BASES},
        alignment_rate => ($qc_metric_result->{PAIR}->{PF_READS_ALIGNED} / $qc_metric_result->{PAIR}->{TOTAL_READS}),
        duplication_rate => $qc_metric_result->{DUPLICATION_RATE},
        first_of_pair_mismatch_rate => $qc_metric_result->{FIRST_OF_PAIR}->{PF_MISMATCH_RATE},
        second_of_pair_mismatch_rate => $qc_metric_result->{SECOND_OF_PAIR}->{PF_MISMATCH_RATE},
    );
    return \%data;
}

#These results only exist from running SomaticValidation models on the same instrument data
sub _add_microarray_vbid_result_metrics_to_hash_ref {
    my $self = shift;
    my $qc_metric_result = shift;
    my $data = shift;

    my $build = Genome::Model::Build->get($qc_metric_result->{build_id});

    my ($vbid) = grep { $_->class eq 'Genome::InstrumentData::VerifyBamIdResult'} $build->results;
    if ($vbid) {
        $data->{'microarray_freemix'} = $vbid->freemix;
        $data->{'microarray_chipmix'} = $vbid->chipmix;
    } else {
        $data->{'microarray_chipmix'} = undef;
        $data->{'microarray_freemix'} = undef;
    }
    return 1;
}

sub _write_metrics_hash_refs_to_tsv_file {
    my $self = shift;
    my $headers = shift;
    my $data = shift;

    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        headers => $headers,
        separator => "\t",
        in_place_of_null_value => 'NA',
        output => $self->output_file,
    );
    unless ($writer) {
        die('Failed to open output file: '. $self->output_file);
    }
    for my $metrics_hash_ref (@$data) {
        $writer->write_one($metrics_hash_ref);
    }
    $writer->output->close;
    return 1;
}

1;
