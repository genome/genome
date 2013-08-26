package Genome::Model::Tools::Dindel::Base;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Dindel::Base {
    is => 'Command',
    is_abstract => 1,
};

sub dindel_install_dir {
    my $self = shift;
    return "/gscmnt/gc2146/info/medseq/dindel";
}

sub dindel_executable {
    my $self = shift;
    return File::Spec->join($self->dindel_install_dir, "binaries", "dindel-1.01-linux-64bit");
}

1;
