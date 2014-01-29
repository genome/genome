package Genome::Model::Tools::Soap::DeNovoAssemble;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Soap::DeNovoAssemble {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
	version => {
	    is => 'Text',
	    doc => 'Version of Soap DeNovo to use',
	    valid_values => ['1.03', '1.04', '1.05', '2.23'],
	},
	config_file => {
	    is => 'Text',
	    doc => 'Config file with run params and input fasta file locations',
	},
	output_dir_and_file_prefix => {
	    is => 'Text',
	    doc => 'Path and common prefix name for output files',
	},
    ],
    has_optional => [
	kmer_size => { #-K
	    is => 'Number',
	    doc => 'K-mer size',
	    default_value => 23,
	    valid_values => [13, 15, 17, 19, 21, 23, 25, 27, 29, 31],
	},
	resolve_repeats => { #-R
	    is => 'Boolean',
	    doc => 'Resolve repeats by reads',
	},
	kmer_frequency_cutoff => { #-d
	    is => 'Number',
	    doc => 'Delete kmers with this frequency or less',
	    valid_values => [1 .. 3], #any value > 0 ??
	},
	cpus => { #-p
	    is => 'Number',
	    doc => 'Number of CPUs to use',
	    default_value => 8,
	    valid_values => [1..8], #any value > 0?
	},
	merge_level => { #-M
	    is => 'Number',
	    doc => 'Strength of merging similar sequences',
	    valid_values => [1..3],
	},
	edge_coverage_cutoff => { #-D
	    is => 'Number',
	    doc => 'Delete contig end sequences with this coverage or less',
	    valid_values => [1, 2], #any value > 0??
	},
	min_scaffold_contig_length => { #-L
	    is => 'Number',
	    doc => 'Exclude contigs less than this length from scaffolding',
	},
    fill_gaps_in_scaffold => { #-F
	    is => 'Boolean',
	    doc => 'Fill gaps in scaffold',
    },
	#TODO .. more options to add
	#-G gapLenDiff(default 50): allowed length difference between estimated and filled gap
	#-u (optional): un-mask contigs with high coverage before scaffolding (default mask)
    ],
};

sub help_brief {
    'Tool to run soap denovo assembler';
}

sub help_detail {
    return <<"EOS"
gmt soap de-novo-assemble --version 1.03 --config-file /gscmnt/111/soap_donovo_assembly/config.txt --output-dir-and-file-prefix /gscmnt/111/soap_donovo_assembly/61BKE
EOS
}

sub execute {
    my $self = shift;

    #TODO - must run on 64 bit

    #validate config file and file contents
    unless ($self->_validate_config_file) {
	$self->error_message("Failed to validate config file");
	return;
    }

    #check output_and_prefix .. make sure output dir exists
    my $output_dir = File::Basename::dirname($self->output_dir_and_file_prefix);
    unless (-d $output_dir) {
	$self->error_message("Invalid output directory: $output_dir does not exist");
	return;
    }

    #eg $ENV{GENOME_SW}/soap/SOAPdenovo-1.04/SOAPdenovo all -s 61BKE_untrimmed.config -K 31 -R -d 1 -o 61BKE_Untrimmed -p 8

    #required params
    my $cmd = $self->path_for_soap_denovo_version.' all -s '.$self->config_file.' -o '.$self->output_dir_and_file_prefix;
    #optional with default values
    $cmd .= ' -K '.$self->kmer_size if $self->kmer_size;
    $cmd .= ' -p '.$self->cpus if $self->cpus;
    $cmd .= ' -M '.$self->merge_level if defined $self->merge_level;
    $cmd .= ' -R '.$self->resolve_repeats if $self->resolve_repeats;
    $cmd .= ' -d '.$self->kmer_frequency_cutoff if $self->kmer_frequency_cutoff;
    $cmd .= ' -D '.$self->edge_coverage_cutoff if $self->edge_coverage_cutoff;
    $cmd .= ' -L '.$self->min_scaffold_contig_length if $self->min_scaffold_contig_length;
    $cmd .= ' -F ' if $self->fill_gaps_in_scaffold;

    $self->debug_message("Running SOAPdenovo with command: $cmd");

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message("Failed to run soap denovo: $@");
        return;
    }

    return 1;
}

sub _validate_config_file {
    my $self = shift;

    unless (-s $self->config_file) {
	$self->error_message("Failed to find config file: ".$self->config_file);
	return;
    }

    #validate input files in config file
    my @missing_inputs;
    my $fh = Genome::Sys->open_file_for_reading($self->config_file);
    while (my $line = $fh->getline) {
	chomp $line;
	if ($line =~ /^[fq]\d+\=/ or $line =~ /^[fq]\=/) {
	    my ($input_file) = $line =~ /\=(\S+)$/;
	    unless (-s $input_file) {
		push @missing_inputs, $input_file;
	    }
	}
    }
    $fh->close;

    if (@missing_inputs) {
	my $message = "Failed to find following input file(s) specified in config file:\n";
	foreach (@missing_inputs) {
	    $message .= "\t".$_."\n";
	}
	$self->error_message("$message");
	return;
    }

    return 1;
}

1;
