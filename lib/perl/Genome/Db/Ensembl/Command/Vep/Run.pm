package Genome::Db::Ensembl::Command::Vep::Run;

use strict;
use warnings;
use Genome;

class Genome::Db::Ensembl::Command::Vep::Run {
    is => 'Genome::Db::Ensembl::Command::Vep::Base',
    has => [
        ensembl_version => {
            is => 'String',
            doc => 'Version of ensembl database to use',
        },
    ],
};

1;

