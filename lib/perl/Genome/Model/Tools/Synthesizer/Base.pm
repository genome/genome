package Genome::Model::Tools::Synthesizer::Base;

use strict;
use warnings;

use Genome;


class Genome::Model::Tools::Synthesizer::Base {
	is          => 'Command::V2',
	is_abstract => 1,
	
};
sub help_detail {
    "Toolkit for downstream analysis of sncRNA sequence data. These commands are setup to run on perl5.12.1";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    unless ($] > 5.010) {
        die 'Bio::DB::Sam requires perl 5.10!';
    }
    require Bio::DB::Sam;
    return $self;
}


1;

