package Genome::Db::Cosmic::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::Cosmic::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'cosmic' 
                },
    ],
    doc => 'install some version of the COSMIC data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db cosmic install 65.1
EOS
}

1;
