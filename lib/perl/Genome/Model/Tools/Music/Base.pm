package Genome::Model::Tools::Music::Base;

use strict;
use warnings;
use Genome;

our $VERSION = $Genome::Model::Tools::Music::VERSION;

class Genome::Model::Tools::Music::Base {
    is => ['Command::V2'],
    is_abstract => 1,
    attributes_have => [
        file_format => {
            is => 'Text',
            is_optional => 1,
        }
    ],
    doc => "Mutational Significance In Cancer (Cancer Mutation Analysis)"
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
    return <<EOS
 Nathan D. Dees, Ph.D.
 Cyriac Kandoth, Ph.D.
 Dan Koboldt, M.S.
 William Schierding, M.S.
 Michael Wendl, Ph.D.
 Qunyuan Zhang, Ph.D.
 Thomas B. Mooney, M.S.           
EOS
}

=cut
sub _doc_credits {
    return ('','None at this time.');
}
=cut

sub _doc_see_also {
    return <<EOS
B<genome-music>(1), B<genome>(1)
EOS
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
}

sub help_detail {
    return "This tool is part of the MuSiC suite. See:\n" 
        . join("\n",shift->_doc_see_also)
        . "\n";
}

1;
