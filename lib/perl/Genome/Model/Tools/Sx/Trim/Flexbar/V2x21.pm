package Genome::Model::Tools::Sx::Trim::Flexbar::V2x21;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Trim::Flexbar::V2x21 {
    is => 'Genome::Model::Tools::Sx::Trim::Flexbar::Base',
    has => { __PACKAGE__->_cmd_properties, },
};

sub executable_path { '/usr/bin/flexbar229'; }
sub version { '2.21'; }

sub input_param_name { 'source'; }
sub input_format_param_value { 'fastq-sanger'; }

sub help_brief { 'Flexbar version 2.21 [commit 229: 2.21 + bug fixes]'; }

our @trim_end_modes = (qw/ ANY LEFT RIGHT LEFT_TAIL RIGHT_TAIL /);
sub _cmd_properties {
    return (
        threads => {
            is => 'Number',
            is_optional => 1,
            default_value => 1,
            doc => 'Number of threads, default: 1',
        },
        barcodes => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file with barcodes, specify (br) to use seperate barcode reads',
        },
        barcode_reads => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta or fastq file with barcode reads, if barcodes not within reads',
        },
        barcode_trim_end => {
            is => 'Text',
            is_optional => 1,
            valid_values => \@trim_end_modes,
            doc => 'Type of barcoding within source reads, see section trim-end types, default: ANY',

        },
        barcode_threshold => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed mismatches and indels per 10 bases for barcode, default: 1.0',
        },
        barcode_min_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap for barcodes (default is length of first barcode)',
        },
        barcode_match => {
            is => 'Number',
            is_optional => 1,
            doc => 'Match score, default: 1',
        },
        barcode_mismatch => {
            is => 'Number',
            is_optional => 1,
            doc => 'Mismatch score, default: -1',
        },
        barcode_gap_cost => {
            is => 'Number',
            is_optional => 1,
            doc => 'Gap score, default: -7',
        },
        barcode_remove => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Remove barcodes within reads based on barcoding parameters',
        },
        adapters => {
            is => 'Text',
            is_optional => 1,
            doc => 'Fasta file with adapters, or barcodes to remove within reads',
        },
        adapter_seq => {
            is => 'Text',
            is_optional => 1,
            doc => 'Single adapter sequence as alternative to adapters option',
        },
        adapter_trim_end => {
            is => 'Text',
            is_optional => 1,
            valid_values => \@trim_end_modes,
            doc => 'Type of alignment for removal, see section trim-end types, default: RIGHT',
        },
        adapter_threshold => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed mismatches and indels per 10 bases for adapter, default: 3.0',
        },
        adapter_min_overlap => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum overlap of adapter and read in base pairs, default: 8',
        },
        adapter_match => {
            is => 'Number',
            is_optional => 1,
            doc => 'Match score, default: 1',
        },
        adapter_mismatch => {
            is => 'Number',
            is_optional => 1,
            doc => 'Mismatch score, default: -1',
        },
        adapter_gap_cost => {
            is => 'Number',
            is_optional => 1,
            doc => 'Gap score, default: -7',
        },
        adapter_no_adapt => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Do not treat parameter min-overlap as adaptive measure, see doc',
        },
        min_readlength => {
            is => 'Number',
            is_optional => 1,
            doc => 'Minimum read length to remain after removal, default: 18',
        },
        max_uncalled => {
            is => 'Number',
            is_optional => 1,
            doc => 'Allowed uncalled bases (N or .) in reads, default: 0',
        },
        pre_trim_front => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim specified number of bases on 5\' end of reads before removal',
        },
        pre_trim_back => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim specified number of bases on 3\' end of reads before removal',
        },
        pre_trim_phred => {
            is => 'Number',
            is_optional => 1,
            doc => 'Trim reads from 3\' end until specified or higher quality reached',
        },
        no_length_dist => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'Prevent length distribution for each read output file',
        },
        removal_tag => {
            is =>'Boolean',
            is_optional => 1,
            doc => 'Tag reads for which adapter or barcode is removed',
        },
        log_level => {
            is => 'Text',
            is_optional => 1,
            doc => 'Print alignments for all or modified reads. One of ALL, MOD, and TAB',
        },
    );
}

1;

