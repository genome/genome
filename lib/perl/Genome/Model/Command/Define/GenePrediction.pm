package Genome::Model::Command::Define::GenePrediction;

use strict;
use warnings;

use Genome;

class Genome::Model::Command::Define::GenePrediction {
    is => 'Command::SubCommandFactory',
    is_abstract => 1,
    doc => 'define a new gene prediction model',
};

# All subclasses of this class have a model define command generated
sub _sub_commands_from { 'Genome::Model::GenePrediction' };
sub _target_base_class { 'Genome::Model::GenePrediction' };

# All generated subcommands inherit from this class
sub _sub_commands_inherit_from { 'Genome::Model::Command::Define::GenePrediction::Helper' };

1;

