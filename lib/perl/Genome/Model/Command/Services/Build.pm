
package Genome::Model::Command::Services::Build;

use strict;
use warnings;

use Genome;
use Command; 

class Genome::Model::Command::Services::Build {
    is => 'Command',
    doc => "build services to be run out of cron or lsf"
};

sub sub_command_sort_position { 98 }

1;

