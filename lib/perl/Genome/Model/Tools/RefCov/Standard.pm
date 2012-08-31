package Genome::Model::Tools::RefCov::Standard;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::RefCov::Standard {
    is => ['Genome::Model::Tools::RefCov'],
};

sub help_brief {
    "The standard mode makes all optional ref-cov outputs available to the user, but none are set by default.  The standard mode can be used to run ref-cov in a customized fashion.",
}

sub execute {
    my $self = shift;
    unless ($self->print_roi_coverage) {
        die('Failed to print ROI coverage!');
    }
    return 1;
}


1;


