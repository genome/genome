package Genome::Model::Tools::Bacterial;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Bacterial {
    is => ['Command'],
    english_name => 'bap',

};

sub help_brief {
    "tools to work for the bacterial annotation pipeline"
}


sub sub_command_sort_position { 16 }

sub help_detail {
    return <<EOS
need to fill out the bacterial help detail
EOS
}

sub command_name {
    'bap' 
}




1;

