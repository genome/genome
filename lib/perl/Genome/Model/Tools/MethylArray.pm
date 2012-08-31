package Genome::Model::Tools::MethylArray;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MethylArray {
    is  => 'Command::Tree',
    doc => 'Toolkit for processing 450k Infinium Methylation Array data',
};


sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt methyl-array ...    
EOS
}


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    return $self;
}



1;