package Genome::Model::Tools::Kmer::Suffixerator;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Kmer::Suffixerator {
    is => 'Genome::Model::Tools::Kmer',
    has => [
        fasta_files => { is => 'Text',
                         doc => 'All the sequence fasta files to perform kmer analysis',
                         is_many => 1
                     },
        index_name => {
            is => 'Text',
            doc => 'The path to the suffix array index',
        },
    ],
    has_optional => [
        parts => { is => 'Integer',
                   doc => 'The number of parts to split suffix in memory(default_value=4).',
                   default_value => 4
               },
        log_file => { is => 'Text', },
    ],
};

sub execute {
    my $self = shift;
    my $genometools_path =$self->genometools_path;
    my @fasta_files = $self->fasta_files;
    my $fasta_files_string = join(' ',@fasta_files);
    #Something went wrong with the 64-bit install or the build of 64-bit binaries
    my $db_flag = '-db';
    if (Genome::Config->arch_os =~ /64/) {
        $db_flag = '-forward-delete-char';
    }
    my $cmd = $genometools_path .' suffixerator -dna -pl -tis -suf -lcp -v -parts '. $self->parts .' '. $db_flag  .' '. $fasta_files_string .' -indexname '. $self->index_name;
    if ($self->log_file) {
        $cmd .= ' > '. $self->log_file;
    };
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => \@fasta_files,
    );
    return 1;
}
