package Genome::Model::Tools::Fasta::Deduplicator;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::Fasta::Deduplicator
{
    is => 'Genome::Model::Tools::Fasta',
    has_input => [
            deduplicated_file   => {
                                    doc => 'file to write contaminations to',
                                    is => 'String',
                                    is_output => 1,
                                    is_optional => 1,
                                },
         ],
};

sub help_brief 
{
    "remove duplicates from a file of reads";
}

sub help_synopsis 
{
    return <<"EOS"
EOS
}

sub create 
{
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return $self;
}

sub execute 
{
    my $self = shift;
    my $fasta_file = $self->fasta_file;
    my $deduplicated_file = ($self->deduplicated_file ? $self->deduplicated_file : $fasta_file . "dedup");
    my $fa_in_io = $self->get_fasta_reader($fasta_file);
    my $fa_out_io = $self->get_fasta_writer($deduplicated_file);
    my %reads;
    
    while (my $seq = $fa_in_io->next_seq) 
    {
        my $val = $seq->id;
        ($fa_out_io->write_seq($seq) and $reads{$val}++) unless ($reads{$val});
    }   
    $self->deduplicated_file($deduplicated_file);

    return 1;
}

1;
