package Genome::Db::CancerGeneList::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::CancerGeneList::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'cancer-gene-list' 
                },
    ],
    doc => 'install some version of the CancerGeneList data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db cancer-gene-list install --species=human --branch=human-1
EOS
}

1;
