package Genome::Model::CwlPipeline::Command;

use Genome;

use strict;
use warnings;

class Genome::Model::CwlPipeline::Command {
    is => 'Command::Tree',
    doc => 'operate on cwl-pipeline models',
};

sub sub_command_category { 'type specific' }

1;

