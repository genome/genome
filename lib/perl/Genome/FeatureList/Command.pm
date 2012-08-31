package Genome::FeatureList::Command;

use strict;
use warnings;

use Genome;


class Genome::FeatureList::Command {
    is => 'Command::Tree',
    has => [
        feature_list => {
            is => 'Genome::FeatureList',
            shell_args_position => 1,
        },
    ],
};

sub help_brief {
    "work with feature-lists"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt feature-list ...    
EOS
}

sub help_detail {                           
    return <<EOS 
A collection of commands to interact with feature-lists.
EOS
}

1;
