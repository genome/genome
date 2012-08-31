package Genome::Model::MutationalSignificance::Command; 

use strict;
use warnings;

use Genome;

class Genome::Model::MutationalSignificance::Command {
    is => 'Command::Tree',
};

#< Help >#
sub help_brief {
    return 'operate on mutational signifigance models/builds';
}

sub help_detail {
    return help_brief();
}
#<>#

#< Command API >#
sub sub_command_category { 'type specific' }

sub _command_name_brief {
    return 'mutational-signifigance';
}
#<>#

1;

