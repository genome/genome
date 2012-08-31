package Genome::Model::Tools::SmrtAnalysis;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis {
    is  => 'Command::Tree',
    doc => 'The Pacific Biosciences secondary analysis software suite.',
};

sub help_brief {
    "Tools to run the SMRT Analysis package from Pacific Biosciences.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt smrt-analysis ...    
EOS
}

sub help_detail {
    return <<EOS 
More information about the SMRT Analysis suite of tools can be found at http://pacbio/ (INTERNAL) or http://pacbio.force.com/DevNet (EXTERNAL).
EOS
}

1;

