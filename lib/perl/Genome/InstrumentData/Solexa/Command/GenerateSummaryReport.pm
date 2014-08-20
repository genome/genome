package Genome::InstrumentData::Solexa::Command::GenerateSummaryReport; #TODO: name this for real

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Solexa::Command::GenerateSummaryReport {
    is => 'Command::V2',
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData::Solexa',
            shell_args_position => 1,
            doc => '', #TODO: write me
        },
    ],
    doc => '', #TODO: write me
};

sub execute {
    my $self = shift;
    my $instrument_data = $self->instrument_data;
    my @alignment_results = Genome::InstrumentData::AlignmentResult->get(instrument_data_id => $instrument_data->id);
    my $counter = 1;
    for my $result (@alignment_results){
        my $instrument_data = $result->instrument_data;
        my %metric_hash = create_metric_hash($result);
        my %attributes_hash = create_attributes_hash($instrument_data);

        print join("\t", $counter,
                         $instrument_data->flow_cell_id,
                         $instrument_data->lane,
                         $attributes_hash{'fwd_run_type'},
                         $attributes_hash{'fwd_read_length'},
                         '', #library (empty on report)
                         $instrument_data->sample->individual_common_name . ' ' . $instrument_data->sample->common_name,
                         '', #library type (empty on report)
                         '', #Loaded [Library] pM
                         $instrument_data->median_insert_size,
                         $metric_hash{'median_insert_size'},
                         $instrument_data->clusters,
                         $metric_hash{'paired_end_base_count'},
                         $instrument_data->fwd_filt_error_rate_avg,
                         $instrument_data->read_1_pct_mismatch,
                         $attributes_hash{'filt_error_rate_avg'},
                         $attributes_hash{'fwd_filt_error_rate_avg'},
                         $attributes_hash{'fwd_filt_aligned_clusters'},
                         $metric_hash{'read_1_pct_aligned'},
                         $attributes_hash{'rev_filt_aligned_clustered_pct'},
                         $metric_hash{'read_2_pct_aligned'},
                         _calculate_avg_qscore(\%attributes_hash, 'fwd', $instrument_data->clusters), #Avg QScore (R1)
                         _calculate_avg_qscore(\%attributes_hash, 'rev', $instrument_data->clusters), #Avg QScore (R2)
                         _calculate_percent_gt_q30(\%attributes_hash, 'fwd', $instrument_data->clusters), # % >Q30 (R1)
                         _calculate_percent_gt_q30(\%attributes_hash, 'rev', $instrument_data->clusters), # % >Q30 (R2)
                         '', #Percent Unknown PF Reads (empty on report)
                         '', #Projects
                         '', #Work Order (id)
                         '', #SRA experiment Accession (empty on report)
                         '', #SRA Run Accession (empty on report
                         '', #SRA Run Status (empty on report)
                  ), "\n";
        $counter++;
    }
    return 1;
}

sub _calculate_avg_qscore {
    my ($attributes_hash, $prefix, $clusters) = @_;
    my $base_quality_sum = $attributes_hash->{$prefix . '_base_quality_sum'};
    my $read_length = $attributes_hash->{'fwd_read_length'};
    return (($base_quality_sum / ($clusters * $read_length)) / 100);
}

sub _calculate_percent_gt_q30 {
    my ($attributes_hash, $prefix, $clusters) = @_;
    my $q30_base_count = $attributes_hash->{$prefix . '_q30_base_count'};
    $DB::single = 1 if $q30_base_count == 0;
    my $read_length = $attributes_hash->{'fwd_read_length'};
    return ($q30_base_count / ($clusters * $read_length));
}

sub create_metric_hash {
    my $r = shift;
    my @metrics = $r->metrics;
    my %metric_hash;
    for my $m (@metrics){
        $metric_hash{$m->metric_name} = $m->metric_value;
    }
    return %metric_hash;
}

sub create_attributes_hash {
    my $instrument_data = shift;
    my @attributes = $instrument_data->attributes;;
    my %attributes_hash;
    for my $a (@attributes){
        $attributes_hash{$a->attribute_label} = $a->attribute_value;
    }
    return %attributes_hash;
}

1;
