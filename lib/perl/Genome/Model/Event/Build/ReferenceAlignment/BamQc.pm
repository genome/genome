package Genome::Model::Event::Build::ReferenceAlignment::BamQc;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::ReferenceAlignment::BamQc {
    is  => ['Genome::Model::Event'],
    has_transient_optional => [
        _alignment_result => { is => 'Genome::SoftwareResult'}
    ],
    doc => 'runs BamQc on the bam(s) produced in the alignment step',
};


sub shortcut {
    my $self = shift;

    my %params = $self->params_for_result;
    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_with_lock(%params);

    if ($result) {
        $self->status_message('Using existing result ' . $result->__display_name__);
        return $self->link_result($result);  #add user and symlink dir to alignment result output dir
    } 
    else {
        return;
    }
}

sub execute {
    my $self = shift;
    my $build = $self->build;

    my %params = (
        $self->params_for_result,
        log_directory => $build->log_directory,
    );

    my $result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->get_or_create(%params);
    $self->link_result($result);

    $self->status_message('Set several metrics to alignment result');
    
    my %qc_metrics;
    map{$qc_metrics{$_->metric_name} = $_->metric_value}$result->metrics;

    my $instr_data = $self->instrument_data;
    my %metrics_to_add;

    if ($instr_data->is_paired_end) {
        my %convert = (
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-FIRST_OF_PAIR',  'PCT_PF_READS_ALIGNED') => 'read_1_pct_aligned',
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-SECOND_OF_PAIR', 'PCT_PF_READS_ALIGNED') => 'read_2_pct_aligned',
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-FIRST_OF_PAIR',  'PF_MISMATCH_RATE')     => 'read_1_pct_mismatch',
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-SECOND_OF_PAIR', 'PF_MISMATCH_RATE')     => 'read_2_pct_mismatch',
        );

        %metrics_to_add = $self->convert_metrics(\%convert, \%qc_metrics, 1);
    }
    else { #Only one set of alignment metrics for single_end instrument data
        my %convert = (
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-UNPAIRED', 'PCT_PF_READS_ALIGNED') => 'read_1_pct_aligned',
            _resolve_metric_name('AlignmentSummaryMetrics', 'CATEGORY-UNPAIRED', 'PF_MISMATCH_RATE')     => 'read_1_pct_mismatch',
        );

        %metrics_to_add = $self->convert_metrics(\%convert, \%qc_metrics, 1);
    }

    my %convert = (
        _resolve_metric_name('InsertSizeMetrics', 'PAIR_ORIENTATION-FR', 'MEDIAN_INSERT_SIZE') => 'median_insert_size',
        _resolve_metric_name('InsertSizeMetrics', 'PAIR_ORIENTATION-FR', 'STANDARD_DEVIATION') => 'sd_insert_size',
    );
    %metrics_to_add = (%metrics_to_add, $self->convert_metrics(\%convert, \%qc_metrics, 0));

    my $align_result = $self->_alignment_result;

    my %align_metrics;
    map{$align_metrics{$_->metric_name} = $_->metric_value}$align_result->metrics;

    for my $metric_key (sort keys %metrics_to_add) {
        if (exists $align_metrics{$metric_key}) {
            $self->warning_message("metric: $metric_key already exist for ".$align_result->id.'. Skip.');
            next;
        }
        $align_result->add_metric(metric_name => $metric_key, metric_value => $metrics_to_add{$metric_key});
    }

    return 1;
}

sub params_for_result {
    my $self  = shift;
    my $build = $self->build;
    my $pp    = $build->processing_profile;

    unless ($self->_alignment_result) {
        my $instrument_data_input = $self->instrument_data_input;
        my ($align_result) = $build->alignment_results_for_instrument_data($instrument_data_input->value);

        unless ($align_result) {
            die $self->error_message('No alignment result found for build: '. $build->id);
        }
        $self->_alignment_result($align_result);
    }

    my $picard_version = $pp->picard_version;

    if ($picard_version < 1.40) {
        my $pp_picard_version = $picard_version;
        $picard_version = Genome::Model::Tools::Picard->default_picard_version;
        $self->warning_message('Given picard version: '.$pp_picard_version.' not compatible to CollectMultipleMetrics. Use default: '.$picard_version);
    }

    my $instr_data  = $self->instrument_data;
    my $er_pileup   = $pp->read_aligner_name =~ /^bwa$/i ? 0 : 1;
    my $read_length = $instr_data->sequencing_platform =~ /^solexa$/i ? 0 : 1;

    return (
        alignment_result_id => $self->_alignment_result->id,
        picard_version      => $picard_version,
        samtools_version    => $pp->samtools_version,
        fastqc_version      => '0.10.0',
        samstat_version     => Genome::Model::Tools::SamStat::Base->default_samstat_version,
        error_rate          => 1,
        error_rate_pileup   => $er_pileup,
        read_length         => $read_length,
        test_name           => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
}

sub link_result {
    my ($self, $result) = @_;

    my $build = $self->build;
    $result->add_user(label => 'uses', user => $build);

    my $align_result = $self->_alignment_result;
    my $link = join('/', $align_result->output_dir, 'bam-qc-'.$result->id);

    if (-l $link) {
        $self->status_message("Already found a symlink here.");
    }
    else {
        Genome::Sys->create_symlink($result->output_dir, $link);
        $align_result->_reallocate_disk_allocation;
    }

    my @users = $align_result->users;
    unless (grep{$_->user eq $result}@users) {
        $align_result->add_user(label => 'uses', user => $result);
    }

    return 1;
}


sub convert_metrics {
    my ($self, $convert, $metrics, $flag) = @_;
    my %convert_metrics;

    for my $metric_name (sort keys %$convert) {
        my $metric_value = $metrics->{$metric_name};
        if ($metric_value) {
            $metric_value = sprintf("%.2f", $metric_value * 100) if $flag;
        }
        else {
            $self->warning_message("Failed to get metric: $metric_name from BamQc output");
            $metric_value = 0;
        }
        $convert_metrics{$convert->{$metric_name}} = $metric_value;
    }

    return %convert_metrics;
}

#convert to BamQc software result metric name
sub _resolve_metric_name {
    my ($type, $key, $label) = @_;
    return sprintf('bam_qc-%s-%s-%s', $type, $key, $label);
}

1;
