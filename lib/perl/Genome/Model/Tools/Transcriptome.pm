package Genome::Model::Tools::Transcriptome;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Transcriptome {
    is  => 'Command::Tree',
    doc => 'Toolkit for transcriptome sequence analysis.',
};

sub help_brief {
    "Assess RNA-Seq, cDNA-Capture, RNA-Cap etc. sequence data against the transcriptome.",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt transcriptome ...    
EOS
}
sub help_detail {
    "The commands in this tree are intended to analyze transcriptome sequence data.  This includes RNA-Seq, RNA-Cap, cDNA-Capture and other library preps.  Analysis would include transcript characterization, expression, splicing, fusion, etc.";
}


1;
