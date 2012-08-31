package Genome::Model::Tools::WuBlast::Blastn;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::WuBlast::Blastn {
    is => 'Genome::Model::Tools::WuBlast',
    has_optional => [
    N => {
        is => 'Integer',
        doc => 'Mismatch score',
    },
    ],
};
sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    # Set output TODO move to class def and use calculate
    $self->output_file( $self->database.'.blast') unless defined $self->output_file;
    unlink $self->output_file if -e $self->output_file;
    return $self;
}

sub _additional_blast_params {
    return (qw/ N /);
}

1;

#$HeadURL$
#$Id$
