package Genome::Model::Event::Build::AmpliconAssembly::Orient;

use strict;
use warnings;

use Genome;

class Genome::Model::Event::Build::AmpliconAssembly::Orient {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    my $orient = Genome::Model::Tools::AmpliconAssembly::Orient->create(
        directory => $self->build->data_directory,
    )
        or return;
    $orient->execute
        or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
