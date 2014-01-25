package Genome::Db::Tgi::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::Tgi::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'tgi' 
                },
    ],
    doc => 'install some data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db tgi install --subsource=cancer-annotation --species=human --branch=human-build37-20130401
genome db tgi install --subsource=misc-annotation --species=human --branch=human-build37-20130113
EOS
}

1;
