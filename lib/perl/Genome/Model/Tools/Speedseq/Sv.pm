package Genome::Model::Tools::Speedseq::Sv;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Speedseq::Sv {
	is => 'Genome::Model::Tools::Speedseq::Base',
	has_param => [

		full_bam_file => {
			is => 'Text',
			doc => 'full BAM file(s) (comma separated list)',
			example_values => ['input1.bam,input2.bam'],			
			tool_param_name => 'B',
			tool_input_file => 1,
		},
		split_read_bam_file => {
                        is => 'Text',
                        doc => 'split reads BAM file(s) (comma separated list)(Auto Generated if absent)',
                        example_values => ['input1.splitters.bam,input2.splitters.bam'],
                        tool_param_name => 'S',
			is_optional => 1,
                        tool_input_file => 1,
                },
		discordant_read_bam_file => {
                        is => 'Text',
                        doc => 'discordant reads BAM file(s) (comma separated list) (Auto Generated if absent)',
                        example_values => ['input1.bam'],
                        tool_param_name => 'D',
                        is_optional => 1,
                        tool_input_file => 1,
                },
		reference_fasta => {
			is => 'Text',
			doc => 'fasta file (indexed with bwa)',
			example_values => ['reference.fa'],
			tool_input_file => 1,
			tool_param_name => 'R',
		},
        	output_prefix => {
        		is => 'Text',
        		doc => 'output prefix [in1.bam]',
        		is_optional => 1,
        		example_values => ['in1.fq','in.realign'],
        		tool_param_name => 'o',
        	},
		threads => {
        		is => 'Integer',
        		doc => 'threads',
       			default_value => 1,
        		is_optional => 1,
        		tool_param_name => 't',
        	},
		exclude_bed_file => {
                        is => 'Text',
                        doc => 'BED file to exclude',
                        is_optional => 1,
			tool_input_file => 1,
                        tool_param_name => 'x',
                },
		genotype_svtyper => {
			is => 'Boolean',
			doc => 'genotype SV breakends with svtyper',
			is_optional => 1,
			tool_param_name => 'g',
		},
		CNVnator_read_depth => {
                        is => 'Boolean',
                        doc => 'calculate read-depth with CNVnator',
                        is_optional => 1,
                        tool_param_name => 'd',
                },
		annotate_vcf => {
                        is => 'Boolean',
                        doc => 'annotate the vcf with VEP',
                        is_optional => 1,
                        tool_param_name => 'A',
                },
		min_sample_weight => {
                        is => 'Integer',
                        doc => 'minimum sample weight for a call',
                        default_value => 4,
			is_optional => 1,
                        tool_param_name => 'm',
                },
		trim_threshold => {
                        is => 'Float',
                        doc => 'trim threshold',
			default_value => 0,
                        is_optional => 1,
                        tool_param_name => 'r',
                },
		temp_directory => {
            		is => 'Text',
            		doc => 'temp directory',
            		example_values => ['./output_prefix.XXXXXXXXXXXX'],
            		is_optional => 1,
            		tool_param_name => 'T',
        	},
		keep_temp_files => {
                        is => 'Boolean',
                        doc => 'keep temporary files',
                        is_optional => 1,
                        tool_param_name => 'k',
                },
		probability_curves => {
                    is => 'Boolean',
                    doc => 'output LUMPY probability curves in VCF',
                    is_optional => 1,
                    tool_param_name => 'P',
                }

	],
};


sub _tool_subcommand_name {
	return 'sv';
}

sub _resolve_input_files {
	my $self = shift;
	my $full_file_list = $self->full_bam_file; 
	my $split_file_list = $self->split_read_bam_file;
	my $discord_file_list = $self->discordant_read_bam_file;
	my @full_params = split(',',$full_file_list);
	my @split_params = split(',',$split_file_list);
	my @discord_params = split(',',$discord_file_list);

	my @final_list = (@full_params, @split_params, @discord_params);

	$self->input_files(\@final_list);

	return (@final_list);
}


