use Genome;
use strict;
use warnings;

package Genome::Model::SmallRna::Command;

class Genome::Model::SmallRna::Command {
    is => 'Command::Tree',
    doc => 'operate on small-rna models',
};

sub sub_command_category { 'type specific' }

1;

