package Genome::Model::Tools::Manta::Config;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Manta::Config {
    is => 'Genome::Model::Tools::Manta::Base',
    has_input => [
        reference_fasta => {
            is => 'Text',
            doc => 'samtools-indexed reference fasta file',
            tool_arg_name => 'referenceFasta',
            tool_input_file => 1,
        },
        working_directory => {
            is => 'Text',
            doc => 'Run script and run output will be written to this directory',
            tool_arg_name => 'runDir',
        },
    ],
    has_optional_input => [
        config_file => {
            is => 'Text',
            doc => 'provide a configuration file to override defaults in global config file',
            tool_arg_name => 'config',
            tool_input_file => 1,
        },
        bams => {
            is => 'Text',
            doc => 'Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample.',
            is_many => 1,
            tool_arg_name => 'bam',
            tool_input_file => 1,
        },
        normal_bam_file => {
            is => 'Text',
            doc => 'Normal sample BAM or CRAM file.',
            tool_arg_name => 'normalBam',
            tool_input_file => 1,
        },
        tumor_bam_file => {
            is => 'Text',
            doc => 'Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted.',
            tool_arg_name => 'tumorBam',
            tool_input_file => 1,
        },
    ],
    has_optional_param => [
        exome => {
            is => 'Boolean',
            doc => 'Set options for WES input: turn off depth filters',
            tool_arg_name => 'exome',
        },
        rna => {
            is => 'Boolean',
            doc => 'Set options for RNA-Seq input: turn off depth filters and do not treat anomalous reads as SV evidence when the proper-pair bit is set.',
            tool_arg_name => 'rna',
        },
        unstranded_rna =>{
            is => 'Boolean',
            doc => 'Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand',
            tool_arg_name => 'unstrandedRNA',
        },
    ],
};

sub _tool_subcommand_name {
    return 'configManta.py';
}


1;
