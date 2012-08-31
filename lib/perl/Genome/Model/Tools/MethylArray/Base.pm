package Genome::Model::Tools::MethylArray::Base;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::MethylArray::Base {
	is          => 'Command::V2',
	is_abstract => 1,
	
};
sub help_detail {
    "Toolkit for processing Methyl 450k Infinium Array data";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    return $self;
}


1;