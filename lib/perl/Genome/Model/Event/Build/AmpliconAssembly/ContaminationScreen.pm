package Genome::Model::Event::Build::AmpliconAssembly::ContaminationScreen;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Event::Build::AmpliconAssembly::ContaminationScreen {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    my $classify = Genome::Model::Tools::AmpliconAssembly::ContaminationScreen->create(
        directory => $self->build->data_directory,
        database => '/gsc/var/lib/reference/set/2809160070/blastdb/blast',
        remove_contaminants => 0, # just run
    )
        or return;
    $classify->execute
        or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
