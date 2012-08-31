package Genome::Model::MetagenomicShotgun::Command;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun::Command {
    is => 'Command::Tree',
};

#< Help >#
sub help_brief {
    return 'operate on metagenomic shotgun models/builds';
}

sub help_detail {
    return help_brief();
}
#<>#

#< Command API >#
sub sub_command_category { 'type specific' }

sub _command_name_brief {
    return 'metagenomic-shotgun';
}
#<>#

1;

