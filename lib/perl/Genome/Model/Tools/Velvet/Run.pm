package Genome::Model::Tools::Velvet::Run;

use strict;
use warnings;

use Genome;

use Data::Dumper;

class Genome::Model::Tools::Velvet::Run {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
	    file_name => {
		is => 'String',
		doc => 'comma separated input file name',
		is_mutable => 1,
	    },
    ],
    has_optional => [
	    #HASH SPECIFIC OPTIONS
	    directory => {
		is          => 'String',
		doc         => 'where assembly should be built (DEFAULT: CURRENT DIRECTORY/velvet_run)',
		default     => 'velvet_run',
		is_mutable  => 1,
	    },
	    hash_length => {
		is          => 'Integer', 
		doc         => 'odd integer of 31 or less (DEFAULT: 27)',
		default     => 31,
	    },
	    file_format => {
		is          => 'String',
		doc         => 'input file format: fasta, fastq, fasta.gz, fastq.gz, eland, gerald. (DEFAULT: fastq)',
		default     => 'fastq',
	    },
	    read_type   => {
		is          => 'String',
		doc         => 'read type: short, shortPaired, short2, shortPaired2, long, longPaired. (DEFAULT: shortPaired)',
		default     => 'shortPaired',
	    },
	    #GRAPH SPECIFIC PARAMS
	    cov_cutoff  => {
		is  => 'Number', 
		doc => 'removal of low coverage nodes AFTER tour bus (DEFAULT: 20)',
		default => 20,
	    },
	    read_trkg   => {
		is  => 'Boolean', 
		doc  => 'tracking of short read positions in assembly (DEFAULT: set for tracking)',
		default => 1,
	    },
	    amos_file   => {
		is  => 'Boolean', 
		doc => 'export assembly to AMOS file (DEFAULT: set for exporting to amos file)',
		default => 1,
	    },
	    exp_cov     => {
		is  => 'Number', 
		doc => 'expected coverage of unique regions (DEFAULT: 98)',
		default => 98,
	    },
	    ins_length  => {
		is  => 'Integer', 
		doc => 'expected distance between two paired end reads (DEFAULT: 260)',
		default => 260,
	    },
	    ins_length2 => {
		is  => 'Integer',
		doc => 'expected distance between two paired-end reads in the second short-read dataset (DEFAULT: no read pairing)',
	    },
	    ins_length_long    => {
		is  => 'Integer', 
		doc => 'expected distance between two long paired-end reads (DEFAULT: no read pairing)',
	    },
	    ins_length_sd      => {
		is  => 'Integer', 
		doc => 'est. standard deviation of respective dataset (DEFAULT: 10% of corresponding length)',
	    },
	    ins_length2_sd     => {
		is  => 'Integer', 
		doc => 'est. standard deviation of respective dataset (DEFAULT: 10% of corresponding length)',
	    },
	    ins_length_long_sd => {
		is  => 'Integer', 
		doc => 'est. standard deviation of respective dataset (DEFAULT: 10% of corresponding length)',
	    },
	    min_contig_lgth       => {
		is  => 'Integer', 
		doc => 'minimum contig length exported to contigs.fa file (DEFAULT: hash length * 2)',
		default => 100,
	    },
	    min_pair_count        => {
		is  => 'Integer',
		doc => 'minimum number of paired end connections to justify the scaffolding of two long contigs (DEFAULT: 10)',
	    },
	    max_branch_length     => {
		is  => 'Integer', 
		doc => 'maximum length in base pair of bubble (DEFAULT: 100)',
	    },
	    max_indel_count       => {
		is  => 'Integer', 
		doc => 'maximum length difference allowed between the two branches of a bubble (DEFAULT: 3)',
	    },
	    max_coverage          => {
		is  => 'Number',
		doc => 'removal of high coverage nodes AFTER tour bus (DEFAULT: no removal)',
	    },
	    max_divergence        => {
		is  => 'Number', 
		doc => 'maximum divergence rate between two branches in a bubble (DEFAULT: 0.2)',
	    },
	    max_gap_count         => {
		is  => 'Integer', 
		doc => 'maximum number of gaps allowed in the alignment of the two branches of a bubble (DEFAULT: 3)',
	    },
	    long_mult_cutoff      => {
		is  => 'Integer',
		doc => 'minimum number of long reads required to merge contigs (DEFAULT: 2)',
	    },
	    out_acefile => {
		is      => 'String', 
		doc     => 'name for output acefile, (DEFAULT: velvet_asm.ace)',
	    },
	    afg_file => {
		is      => 'String',
		doc     => 'afg file name',
	    },
	    fastq_file => {
		is => 'String',
		doc => 'fastq file to build ace',
	    },
	    time        => {
		is      => 'String',
		doc     => 'timestamp inside acefile, must be sync with phd timestamp',
	    }, 

     ],
};

sub help_brief {
    "Tool to run complete velvet assembly",
}

sub help_synopsis {
    return <<"EOS"
gmt velvet run --file-names <input file names>
EOS
}

sub help_detail {
    return <<EOS
Tool to run complete velvet assembly
EOS
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);

    return $self;
}


sub execute {
    my $self = shift;

    unless (-s $self->file_name) {
	$self->error_message("Input file is either missing or empty");
	return;
    }

    my $to_hash = Genome::Model::Tools::Velvet::Hash->create(
							     file_name => $self->file_name,
							     directory => $self->directory,
							     hash_length => $self->hash_length,
							     read_type => $self->read_type,
							     file_format => $self->file_format,
							     version => $self->version,
							     );
    unless ($to_hash->execute) {
	$self->error_message("Failed to execute velvet hash");
	return;
    }

    my $to_graph = Genome::Model::Tools::Velvet::Graph->create(
							       directory => $self->directory,
							       amos_file => $self->amos_file,
							       read_trkg => $self->read_trkg,
							       exp_cov => $self->exp_cov,
							       cov_cutoff => $self->cov_cutoff,
							       ins_length => $self->ins_length,
							       version => $self->version,
							       );
    unless ($to_graph->execute) {
	$self->error_message("Failed to execute velvet graph");
	return;
    }

    if ($self->file_format eq 'fastq') {

	mkdir ($self->directory.'/edit_dir', 0775);

	my $to_ace = Genome::Model::Tools::Velvet::ToAce->create(
								 afg_file => $self->_resolve_afg_file,
								 seq_file => $self->directory.'/Sequences',
								 time => $self->_resolve_time,
								 out_acefile => $self->_resolve_ace_out,
								 );

	unless ($to_ace->execute) {
	    $self->error_message("Failed to execute velvet to ace");
	    return;
	}

	my $contigs_fa_file = $self->directory.'/contigs.fa';
	unless (-s $contigs_fa_file) {
	    $self->error_message("Failed to find contigs.fa file: $contigs_fa_file");
	    return;
	}
	#create velvet std output files
#	my $ec = Genome::Model::Tools::Velvet::CreateStdoutFiles->execute(
	my $ec = Genome::Model::Tools::Velvet::StandardOutputs->execute(
#	    input_fastq_file => $self->file_name,
	    directory => $self->directory,
	    );
	unless ($ec) {
	    $self->error_message("Failed to create submission files");
	    return;
	}
    }
    else {
	$self->error_message("Ace file will not be created since input is not a fastq file");
	return;
    }
    
    return 1;
}

sub _resolve_ace_out {
    my $self = shift;
    unless ($self->out_acefile) {
	return $self->directory.'/edit_dir/velvet_asm.ace';
    }
    return $self->ace_out;
}

sub _resolve_time {
    my $self = shift;
    unless ($self->time) {
	my $time = `date "+%a %b %e %T %Y"`;
	chomp $time;
	return $time;
    }
    return $self->time;
}

sub _resolve_afg_file {
    my $self = shift;
    my $afg_file;
    unless ($afg_file = $self->afg_file) {
	$afg_file = $self->directory.'/velvet_asm.afg';
	unless (-s $afg_file) {
	    $self->error_message("Afg file does not exist: $afg_file");
	    return;
	}
    }
    return $afg_file;
}

1;
