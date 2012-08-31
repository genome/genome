package Genome::Model::Tools::Music::Bmr;
use warnings;
use strict;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION; 

class Genome::Model::Tools::Music::Bmr {
    #is  => ['Command::Tree','Genome::Model::Tools::Music::Base'],
    is => ['Command::Tree'],
    doc => "Calculate gene coverages and background mutation rates."
};

sub _doc_copyright_years {
    (2010,2011);
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
    return " Cyriac Kandoth, Ph.D.";
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
B<genome-music-bmr-calc-covg>(1),
B<genome-music-bmr-calc-bmr>(1),
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
