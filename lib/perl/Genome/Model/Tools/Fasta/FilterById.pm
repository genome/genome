package Genome::Model::Tools::Fasta::FilterById;

use strict;
use warnings;

use Genome;
use Workflow;
use Bio::SeqIO;

class Genome::Model::Tools::Fasta::FilterById
{
    is => 'Genome::Model::Tools::Fasta',
    has_input => [
            filter_list => {
                                    doc => 'list of id\'s to filter from fasta_file',
                                    is => 'String',
                                },
            output_file =>  {
                                doc => 'path of output file',
                                is => 'String',
                                is_output => 1,
                                is_optional => 1,
                            },
         ],
};

sub help_brief 
{
    "removes reads from a fasta file based on list of id's";
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
    my $output_file = ($self->output_file ? $self->output_file : $fasta_file . ".clean");
    print "FILTER ENTERED WITH FASTA $fasta_file, filter list " . $self->filter_list . " AND FILTER $output_file\n";
    print "exists:  " . (-e $output_file ? 'yes' : 'no') . "\n";
    -e $output_file and die("$output_file already exists and cannot be written to");

    my $filter_list = new IO::File $self->filter_list;
    chomp(my @ids = <$filter_list>);
    my %ids = map{$_=>$_} @ids; 

    my $fa_in_io = $self->get_fasta_reader($fasta_file);
    my $fa_out_io = $self->get_fasta_writer($output_file);

    while (my $seq = $fa_in_io->next_seq) 
    { 
        my $val = $seq->id;
        $fa_out_io->write_seq($seq) unless ($ids{$val}); 
    }

    print "FILTER completed\n";
    $self->output_file($output_file);
    return 1;
}

1;
