package Genome::Db::DbVar::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::DbVar::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'dbvar' 
                },
    ],
    doc => 'install some version of the DbVar data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db db-var install --species=human --branch=human-build37
EOS
}

1;
