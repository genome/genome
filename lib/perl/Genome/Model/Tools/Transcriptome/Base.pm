package Genome::Model::Tools::Transcriptome::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Transcriptome::Base {
	is => 'Command::V2',
	is_abstract => 1,
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    unless ($self) { return; }
    unless ($] > 5.010) {
        die 'Bio::DB::Sam requires perl 5.10 or greater!';
    }
    require Bio::DB::Sam;
    return $self;
}

1;
