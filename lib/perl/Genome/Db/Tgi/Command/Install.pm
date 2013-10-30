package Genome::Db::Tgi::Command::Install;
use strict;
use warnings;
use Genome;

class Genome::Db::Tgi::Command::Install {
    is => 'Genome::Db::Command::InstallFromGitRepo',
    has => [
        source  => { is => 'Text',
                    is_constant => 1,
                    value => 'cosmic' 
                },
    ],
    doc => 'install some data set on the local GMS for which there is a branch in the canonical repository'
};

sub help_synopsis {
    return <<EOS
genome db tgi install cancer-annotation/human/build37-20131010.1
genome db tgi install misc-annotation/human/build37-20130113.1
EOS
}

1;
