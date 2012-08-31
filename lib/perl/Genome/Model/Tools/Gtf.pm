package Genome::Model::Tools::Gtf;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf {
    is  => 'Command::Tree',
    doc => 'Tools to work with gtf format annotation files',
};

sub help_brief {
    "Tools to work with gtf format annotation files.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt gtf ...    
EOS
}

sub help_detail {
    return <<EOS 
"More information about the gtf format can be found here: http://genome.ucsc.edu/FAQ/FAQformat.html#format4"
EOS
}

1;

