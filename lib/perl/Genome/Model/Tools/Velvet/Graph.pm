package Genome::Model::Tools::Velvet::Graph;

use strict;
use warnings;

use Genome;
use File::Copy;

my %properties = (
    directory   => {
        is  => 'String', 
        doc => 'working directory name, must contain Roadmaps and Sequences (default: ./velvet_run)',
        default => 'velvet_run',
    },
    cov_cutoff  => {
        is  => 'Number', 
        doc => 'removal of low coverage nodes AFTER tour bus (default: no removal)',
    },
    read_trkg   => {
        is  => 'Boolean', 
        doc  => 'tracking of short read positions in assembly (default: no tracking)',
    },
    amos_file   => {
        is  => 'Boolean', 
        is_input => 1,
        doc => 'export assembly to AMOS file (default: no export)',
    },
    exp_cov     => {
        is  => 'Number', 
        doc => 'expected coverage of unique regions (default: no long read resolution)',
    },
    ins_length  => {
        is  => 'Integer', 
        doc => 'expected distance between two paired end reads (default: no read pairing)',
    },
    ins_length2 => {
        is  => 'Integer', 
        doc => 'expected distance between two paired-end reads in the second short-read dataset (default: no read pairing)',
    },
    ins_length_long    => {
        is  => 'Integer', 
        doc => 'expected distance between two long paired-end reads (default: no read pairing)',
    },
    ins_length_sd      => {
        is  => 'Integer', 
        doc => 'est. standard deviation of respective dataset (default: 10% of corresponding length)',
    },
    ins_length2_sd     => {
        is  => 'Integer', 
        doc => 'est. standard deviation of respective dataset (default: 10% of corresponding length)',
    },
    ins_length_long_sd => {
        is  => 'Integer', 
        doc => 'est. standard deviation of respective dataset (default: 10% of corresponding length)',
    },
    min_contig_lgth       => {
        is  => 'Integer', 
        doc => 'minimum contig length exported to contigs.fa file (default: hash length * 2)',
    },
    min_pair_count        => {
        is  => 'Integer',
        doc => 'minimum number of paired end connections to justify the scaffolding of two long contigs (default: 10)',
    },
    max_branch_length     => {
        is  => 'Integer', 
        doc => 'maximum length in base pair of bubble (default: 100)',
    },
    max_indel_count       => {
        is  => 'Integer', 
        doc => 'maximum length difference allowed between the two branches of a bubble (default: 3)',
    },
    max_coverage          => {
        is  => 'Number',
        doc => 'removal of high coverage nodes AFTER tour bus (default: no removal)',
    },
    max_divergence        => {
        is  => 'Number', 
        doc => 'maximum divergence rate between two branches in a bubble (default: 0.2)',
    },
    max_gap_count         => {
        is  => 'Integer', 
        doc => 'maximum number of gaps allowed in the alignment of the two branches of a bubble (default: 3)',
    },
    long_mult_cutoff      => {
        is  => 'Integer',
        doc => 'minimum number of long reads required to merge contigs (default: 2)',
    },
);


class Genome::Model::Tools::Velvet::Graph {
    is           => 'Genome::Model::Tools::Velvet::Base',
    has_optional => [%properties],
};        


sub help_brief {
    "This tool runs velvetg",
}


sub help_synopsis {
    return <<"EOS"
gmt velvet graph 
EOS
}


sub help_detail {
    return <<EOS
Velvetg is the core of Velvet where the de Bruijn graph is built then manipulated.
The required input files are : Roadmaps and Sequences. The output files are: contigs.fa, 
stats.txt, LastGraph, and velvet_asm.afg (if "-amos_file yes"). 
EOS
}


sub create {
    my $class = shift;
    
    my $self  = $class->SUPER::create(@_);
    my $dir   = $self->directory;

    if (-d $dir) {
        $self->status_message("velvetg will run in directory: $dir");
    }
    else {
        $self->error_message("Need velveth output directory: $dir or Run velveth first");
        return;
    }

    for my $file (map{$dir."/$_"}qw(Roadmaps Sequences)) {
        unless (-s $file) {
            $self->error_message("Need $file to proceed velvetg");
            return;
        }
    }
    
    return $self;
}


sub execute {
    my $self = shift;
    
    my $dir     = $self->directory;
    my $command = $self->resolve_version.' '.$dir;

    for my $property_name (keys %properties) {
        next if $property_name eq 'directory';
        if ($self->can($property_name) and defined $self->$property_name) {
            $command .= " -$property_name ";
            if ($property_name =~ /^(read_trkg|amos_file)$/) {
                $command .= 'yes';
            }
            else {
                $command .= $self->$property_name;
            }
        }
    }
        
    my $stats = $dir.'/stats.txt';
    my $prev  = $stats.'.prev';
    move $stats, $prev if -s $stats;
    
    system $command;
    #Very wierd here, all return value for velvetg run is 256
    
    unless (-s $stats) {
        $self->error_message("velvetg failed.");
        move $prev, $stats if -s $prev;
        return;
    }

    return 1;
}


1;
#$HeadURL$
#$Id$

