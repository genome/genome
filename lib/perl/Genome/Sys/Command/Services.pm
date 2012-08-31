
package Genome::Sys::Command::Services;

use strict;
use warnings;

use Genome;
use Command; 

class Genome::Sys::Command::Services {
    is => 'Command',
    doc => "services to be run out of cron or other task scheduler"
};

sub sub_command_sort_position { 98 }

1;

