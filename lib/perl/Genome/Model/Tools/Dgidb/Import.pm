package Genome::Model::Tools::Dgidb::Import;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Dgidb::Import {
    is => 'Genome::Model::Tools::Dgidb',
    has => [],
};

sub help_brief {
    'Import druggable gene datasource into the database'
}

sub help_synopsis {
    return <<EOS
gmt dgidb import ...
EOS
}

1;
