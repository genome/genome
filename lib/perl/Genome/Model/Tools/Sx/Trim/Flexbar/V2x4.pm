package Genome::Model::Tools::Sx::Trim::Flexbar::V2x4;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Trim::Flexbar::V2x4 {
    is => 'Genome::Model::Tools::Sx::Trim::Flexbar::Base',
    has => { __PACKAGE__->_cmd_properties, },
};

sub executable_path { '/usr/bin/flexbar2.4' }
sub version { '2.4'; }

sub input_param_name { 'reads'; }
sub input_format_param_value { 'sanger'; }

sub help_brief { 'Flexbar version 2.4'; }

our @trim_end_modes = (qw/ ANY LEFT RIGHT LEFT_TAIL RIGHT_TAIL /);
sub _cmd_properties {
    return (
        threads => {
            is => 'Number',
            is_optional => 1,
            default_value => 1,
            doc => 'Number of threads to employ. Default 1.',
        },
        barcodes => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file with barcodes for demultiplexing that may contain N.',
        },
        barcode_reads => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta/q file composed of separate barcode reads for detection.',
        },
        barcode_trim_end => {
            is => 'Text',
            is_optional => 1,
            valid_values => \@trim_end_modes,
            doc => 'Type of detection, see section trim-end modes. Default: ANY.',
        },
        barcode_min_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap of barcode and read. Default: barcode length.',
        },
        barcode_threshold => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed mismatches and gaps per 10 bases overlap. Default: 1.0.',
        },
        barcode_unassigned => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Include unassigned reads in output generation.',
        },
        adapters => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file with adapters for removal that may contain N.',
        },
        adapter_seq => {
            is => 'Text',
            is_optional => 1,
            doc => 'Single adapter sequence as alternative to adapters option.',
        },
        adapter_trim_end => {
            is => 'Text',
            is_optional => 1,
            valid_values => \@trim_end_modes,
            doc => 'Type of removal, see section trim-end modes. Default: RIGHT.',
        },
        adapter_min_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap of adapter and read sequence. Default: 1.',
        },
        adapter_threshold => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed mismatches and gaps per 10 bases overlap. Default: 3.0.',
        },
        max_uncalled => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed uncalled bases (N or .) for each read. Default: 0.',
        },
        pre_trim_left => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim given number of bases on 5\' read end before detection.',
        },
        pre_trim_right => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim specified number of bases on 3\' end prior to detection.',
        },
        pre_trim_phred => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim 3\' end until specified or higher quality reached.',
        },
        post_trim_length => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim to specified read length from 3\' end after removal.',
        },
        min_read_length => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum read length to remain after removal. Default: 18.',
        },
        log_level => {
            is => 'Text',
            is_optional => 1,
            doc => 'Print valid optimal read alignment. One of ALL, MOD, and TAB.',
        },
        length_dist => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Generate length distribution for read output files.',
        },
        removal_tags => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Tag reads that are subject to adapter or barcode removal.',
        },
    );
}

1;

