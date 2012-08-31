package Genome::Model::Tools::Sam::Deduplicator;

use strict;
use warnings;

use Genome;
use Workflow;
use File::Basename;
use IO::File;
use Bio::SeqIO;

class Genome::Model::Tools::Sam::Deduplicator
{
    is => 'Genome::Model::Tools::Sam',
    has_input => [
            sam_file            => {
                                    doc         => 'input file to deduplicate',
                                    is          => 'String',
                                   },
            deduplicated_file   => {
                                    doc         => 'file to write contaminations to',
                                    is          => 'String',
                                    is_output   => 1,
                                    is_optional => 1,
                                   },
            deduplicated_count  => {
                                    doc         => '# of reads deduplicated',
                                    is          => 'Integer',
                                    is_output   => 1,
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
    my $sam_file = $self->sam_file;
    my $deduplicated_file = ($self->deduplicated_file ? $self->deduplicated_file : $sam_file . "dedup");
    unlink ($deduplicated_file) if (-e $deduplicated_file);
    my $sam_fh = Genome::Sys->open_file_for_reading($sam_file) or return;
    my $dd_fh = Genome::Sys->open_file_for_writing($deduplicated_file) or return;
    my $dedup_count = 0;


    my %reads;
    
    while (my $sam = $sam_fh->getline)
    { 
        #my ($chr, $pos, $cns_qual, $snp_qual, $map_qual, $rd_depth) = 
        $sam =~ m/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+/;    
        
        ($dd_fh->print($sam) and $reads{$10}++) unless ($reads{$10} and ++$dedup_count);
    }   

    $self->deduplicated_file($deduplicated_file);
    $self->deduplicated_count($dedup_count);

    return 1;
}

1;
