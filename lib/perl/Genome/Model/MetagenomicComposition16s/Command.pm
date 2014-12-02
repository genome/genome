package Genome::Model::MetagenomicComposition16s::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicComposition16s::Command {
    is => 'Command::Tree',
    doc => 'operate on metagenomic-composition-16s models',
};

sub sub_command_category { 'type specific' }

sub _command_name_brief {
    return 'metagenomic-composition-16s';
}

1;
