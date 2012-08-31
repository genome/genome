package Genome::Model::Tools::Velvet::Hash;

use strict;
use warnings;

use POSIX;
use Genome;

class Genome::Model::Tools::Velvet::Hash {
    is           => 'Genome::Model::Tools::Velvet::Base',
    has => [
         file_name  => {
            is      => 'String', 
            doc     => 'input file name',
        }
    ],
    has_optional => [
        directory   => {
            is      => 'String', 
            doc     => 'directory name for output files, default is ./velvet_run',
            default => 'velvet_run',
        },
        hash_length => {
            is      => 'Integer', 
            doc     => 'odd integer (if even, will be decremented) <= 31, if above, will be reduced, default: 23',
            default => 23,
        },
        file_format => {
            is      => 'String',
            doc     => 'input file format: fasta, fastq, fasta.gz, fastq.gz, eland, gerald. default: fasta',
            default => 'fasta',
        },
        read_type   => {
            is      => 'String',
            doc     => 'read type: short, shortPaired, short2, shortPaired2, long, longPaired. default: short',
            default => 'short',
        },
     ],
};
        

sub help_brief {
    "This tool runs velveth",
}


sub help_synopsis {
    return <<"EOS"
gmt velvet hash --file-name name [--directory dir --hash-length 21 --file-format fastq --read-type short]
EOS
}


sub help_detail {
    return <<EOS
Velveth constructs the dataset for the following program, velvetg, and
indicate to the system what each sequence file represents.Velveth takes 
in a number of sequence files, produces a hashtable, then outputs two files 
in an output directory (creating it if necessary), Sequences and Roadmaps, 
which are necessary to velvetg.
EOS
}


sub create {
    my $class = shift;
    
    my $self  = $class->SUPER::create(@_);
    my $dir   = $self->directory;

    unless (-s $self->file_name) {
	$self->error_message("Input file: ".$self->file_name." not existing or is empty");
	return;
    }
    
    if (-d $dir) {
        $self->warning_message("velveth will overwrite output in directory: $dir");
    }
    else {
        mkdir $dir, 0777;
        unless (-d $dir) {
            $self->error_message("Fail to create output directory: $dir");
            return;
        }
    }
        
    return $self;
}


sub execute {
    my $self = shift;
        
    my $command = sprintf(
        '%s %s %d -%s -%s %s',
        $self->resolve_version,
        $self->directory,
        $self->hash_length,
        $self->file_format,
        $self->read_type,
        $self->file_name,
    );
    
    if (system $command) {
        $self->error_message("$command failed.");
        return;
    }

    return 1;
}


1;
#$HeadURL$
#$Id$

