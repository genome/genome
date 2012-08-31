package Genome::Model::DeNovoAssembly::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::DeNovoAssembly::Command {
    is => 'Command::Tree',
    doc => 'operate on de novo assembly models/builds',
};

sub help_detail {
    return help_brief();
}

sub sub_command_category { 'type specific' }

1;
