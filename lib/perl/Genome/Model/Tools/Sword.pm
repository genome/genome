package Genome::Model::Tools::Sword;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sword {
    is  => 'Command::Tree',
    doc => 'Tools to run Sword on RNA-seq BAM files',
};

sub help_brief {
    "Run Sword on BAM files",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt sword ...    
EOS
}
sub help_detail {
    "These tools will run Sword on RNA-seq BAM files.";
}


1;
