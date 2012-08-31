package Genome::Model::Tools::BamQc;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BamQc {
    is  => 'Command::Tree',
    doc => 'Toolkit for evaluating quality of BAM alignments',
};

sub help_brief {
    "Assess quality of BAM alignments.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt bam-qc ...    
EOS
}
sub help_detail {
    "These commands are intended to provide a robust analysis of BAM alignments combining many individual tools into a comprehensive evaluation.";
}


1;
