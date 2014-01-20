package Genome::Db::Ucsc::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::Ucsc::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'ucsc' 
                },
    ],
    doc => 'install some version of the UCSC data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db ucsc install --species=human --branch=human-build37
EOS
}

1;
