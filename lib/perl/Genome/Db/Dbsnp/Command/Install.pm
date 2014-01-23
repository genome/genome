package Genome::Db::Dbsnp::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::Dbsnp::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'dbsnp' 
                },
    ],
    doc => 'install some version of the COSMIC data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db dbsnp install --species=human --branch=human-build37-132
EOS
}

1;
