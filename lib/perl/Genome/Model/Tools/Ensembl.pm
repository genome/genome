package Genome::Model::Tools::Ensembl;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Ensembl {
    is  => 'Command::Tree',
    doc => 'Analysis tools that rely on the locally installed Ensembl API/MySQL database.',
};

sub help_brief {
    "Tools to work with the local Ensembl API.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt ensembl ...    
EOS
}

sub help_detail {
    return <<EOS 

EOS
}


1;

