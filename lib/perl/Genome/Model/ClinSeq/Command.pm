use Genome;
use strict;
use warnings;

package Genome::Model::ClinSeq::Command;

class Genome::Model::ClinSeq::Command {
    is => 'Command::Tree',
    doc => 'operate on clin-seq models',
};

sub sub_command_category { 'type specific' }

1;

