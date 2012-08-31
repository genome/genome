use Genome;
use strict;
use warnings;

package Genome::Model::ReferenceAlignment::CommandCreateMetrics;

class Genome::Model::ReferenceAlignment::Command::CreateMetrics {
    is => 'Command',
    doc => 'create metrics for reference alignment models',
};

sub sub_command_category { 'type specific' }

1;

