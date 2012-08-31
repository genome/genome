package Genome::Model::Event::Build::AmpliconAssembly::Collate;

use strict;
use warnings;

use Genome;

require Genome::Model::Tools::AmpliconAssembly::Collate;

class Genome::Model::Event::Build::AmpliconAssembly::Collate {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    my $collate = Genome::Model::Tools::AmpliconAssembly::Collate->create(
        directory => $self->build->data_directory,
    )
        or return;
    $collate->execute
        or return;

    return 1;
}

#$HeadURL$
#$Id$
