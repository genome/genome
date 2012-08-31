package Genome::Model::Tools::WuBlast::Blastx;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::WuBlast::Blastx {
    is => 'Genome::Model::Tools::WuBlast',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    my ($basename,$dirname) = File::Basename::fileparse($self->query_file);
    unless ($self->output_directory) {
        $self->output_directory($dirname);
    }
    unless ($self->output_file) {
        $self->output_file( $self->output_directory .'/'. $basename .'.blast');
    }
    unlink $self->output_file if -e $self->output_file;
    return $self;
}

1;

