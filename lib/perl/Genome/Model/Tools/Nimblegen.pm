
package Genome::Model::Tools::Nimblegen;

use strict;
use warnings;

use Genome;
use Command; 

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
);

sub sub_command_sort_position { 1_000_000 }

sub help_brief {
    "Tools for working with Nimblegen Solid Phase Capture Arrays"
}

sub help_detail {
    return <<"EOS"
These commands are for generating and evaluating Nimblegen Capture Arrays used in Validation of detected variants    
EOS
}

1;

