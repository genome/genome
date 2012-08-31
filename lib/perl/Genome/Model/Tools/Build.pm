package Genome::Model::Tools::Build;

use strict;
use warnings;

use Genome;
use File::Basename;

class Genome::Model::Tools::Build {
    is => 'Command',
    has => [
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools that operate on Builds.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt build ...
EOS
}

sub help_detail {                           
    return <<EOS 
Tools that operate on Builds.
EOS
}



1;

