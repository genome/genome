package Genome::Model::Tools::RefSeq;

use strict;
use warnings;

use Genome;            

class Genome::Model::Tools::RefSeq {
    is => 'Command',
    has => [ ],
};

sub sub_command_sort_position { 15 }

sub help_brief {
    'Tools for working with RefSeq Fasta files.'
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools ref-seq ...
EOS
}

sub xhelp_detail {                           
    return <<EOS 
EOS
}

1;

