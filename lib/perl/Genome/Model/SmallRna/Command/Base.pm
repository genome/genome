package Genome::Model::SmallRna::Command::Base;

use strict;
use warnings;

use Genome;


class Genome::Model::SmallRna::Command::Base {
	is          => 'Command::V2',
	is_abstract => 1,
	
};
sub help_detail {
    "These commands require at least perl 5.10.";
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
