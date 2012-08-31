package Genome::Model::Tools::DetectVariants2;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2 {
    is => 'Command::Tree',
};

sub sub_command_sort_position { 12 }

sub help_brief {
    "Tools to detect variants and/or filter their results.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt detect-variants2 ...    
EOS
}

sub help_detail {                           
    return <<EOS 
Tools to detect variants and/or filter their results.
EOS
}

1;

