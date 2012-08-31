package Genome::Model::Tools::SmrtAnalysis::SaxonB;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::SmrtAnalysis::SaxonB {
    is => ['Genome::Model::Tools::SmrtAnalysis::Base'],
    has_input => [
        xsl => { doc => 'Stylesheet file' },
        input => { doc => 'Initial source document' },
        output => { doc => 'Output file or directory', is_output => 1, },
    ],
};


sub execute {
    my $self = shift;

    my $cmd = $self->analysis_bin .'/saxonb9';
    $cmd .= ' -xsl:'. $self->xsl .' -s:'. $self->input .' -o:'. $self->output;
    $self->shellcmd(
        cmd => $cmd,
        input_files => [$self->input],
        output_files => [$self->output],
        skip_if_output_is_present => 0,
    );

    return 1;
}


1;
