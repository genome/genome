package Genome::Model::Tools::EpitopePrediction;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::EpitopePrediction::Base {
	is          => 'Command::V2',
	is_abstract => 1,
	
};
sub help_detail {
    "These commands are setup to run epitope prediction";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    return $self;
}


1;
