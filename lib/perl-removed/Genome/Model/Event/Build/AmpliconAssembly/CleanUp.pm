package Genome::Model::Event::Build::AmpliconAssembly::CleanUp;

use strict;
use warnings;

use Genome;
      
class Genome::Model::Event::Build::AmpliconAssembly::CleanUp {
    is => 'Genome::Model::Event',
};

sub execute {
    my $self = shift;

    my $amplicons = $self->build->get_amplicons
        or return;
    
    for my $amplicon ( @$amplicons ) {
        $amplicon->remove_unneeded_files;
    }

    #print $self->build->data_directory,"\n"; <STDIN>;

    return 1;
}

1;

#$HeadURL$
#$Id$
