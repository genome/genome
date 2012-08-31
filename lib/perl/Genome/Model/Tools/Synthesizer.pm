package Genome::Model::Tools::Synthesizer;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Synthesizer {
    is  => 'Command::Tree',
    doc => 'Synthesizer for analyzing sncRNA sequence data',
};

sub help_brief {
    "Toolkit for downstream analysis of sncRNA sequence data ",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt synthesizer ...    
EOS
}

sub help_detail {
    return <<EOS 
Toolkit for downstream analysis of sncRNA sequence data
EOS
}

1;
