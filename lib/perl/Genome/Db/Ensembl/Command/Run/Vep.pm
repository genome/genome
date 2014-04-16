package Genome::Db::Ensembl::Command::Run::Vep;

use strict;
use warnings;
use Genome;

class Genome::Db::Ensembl::Command::Run::Vep {
    is => 'Genome::Db::Ensembl::Command::Run::Base',
    has => [
        ensembl_version => {
            is => 'String',
            doc => 'Version of ensembl database to use',
        },
    ],
};

1;

