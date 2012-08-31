package Genome::Model::SmallRna::Command::Base;

use strict;
use warnings;

use Genome;


class Genome::Model::SmallRna::Command::Base {
	is          => 'Command::V2',
	is_abstract => 1,
	
};
sub help_detail {
    "These commands are setup to run perl5.10 or /usr/bin/perl";
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    unless ($] > 5.010) {
        die 'Bio::DB::Sam requires perl 5.10! Consider using /usr/bin/perl instead of gmt for all small-rna commands!';
    }
    require Bio::DB::Sam;
    return $self;
}


1;
