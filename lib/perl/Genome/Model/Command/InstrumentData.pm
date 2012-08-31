# FIXME ebelter
#  Long: Remove or update to use inputs as appropriate.
#
package Genome::Model::Command::InstrumentData;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::InstrumentData {
    is => 'Genome::Model::Command::BaseDeprecated',
    doc => "assign, list, dump instrument data related to a model",
};

sub sub_command_sort_position { 6 }

1;
