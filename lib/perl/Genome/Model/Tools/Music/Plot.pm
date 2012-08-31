package Genome::Model::Tools::Music::Plot;
use warnings;
use strict;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION; 

class Genome::Model::Tools::Music::Plot {
    is => ['Command::Tree'],
    doc => "Generate relevant plots and visualizations for MuSiC."
};

sub _doc_copyright_years {
    (2010,2012);
}

sub _doc_license {
    my $self = shift;
    my (@y) = $self->_doc_copyright_years;  
    return <<EOS
Copyright (C) $y[0]-$y[1] Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the 
associated LICENSE file in this distribution.
EOS
}

sub _doc_authors {
    return <<EOS
 Cyriac Kandoth, Ph.D.
 Charles Lu, Ph.D.
 Chris Miller, Ph.D.
 Nathan D. Dees, Ph.D.
EOS
}

=cut
sub _doc_credits {
    # used to compose man pages
    # (return a list of strings)
    return ('','None at this time.');
}
=cut

sub _doc_see_also {
    return <<EOS
B<genome-music>(1),
B<genome>(1)
EOS
}

sub _doc_manual_body {
    return shift->help_detail;
}

sub help_detail {
    return "These tools are part of the MuSiC suite.\n";
}

1;
