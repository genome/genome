package Genome::Model::Tools::Consed::NavToAcenav;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Consed::NavToAcenav {
    is => 'Command',
    has => [
    nav => {
        is => 'Text',
        doc => 'The nav file to convert to an ace nav',
    },
    acenav => {
        is => 'Text',
        doc => 'The new ace nav file to save'
    },
    acefile => {
        is => 'Text',
        doc => 'Acefile full path to add to nav file',
    },
    ],
};

#############

sub help_brief {
    return 'Add acefile name to a nav file';
}

sub help_detail {
    return help_brief();
}

#############

sub execute {
    my $self = shift;

    my $reader = Genome::Model::Tools::Consed::Navigation::Reader->new(
        input => $self->nav,
    )
        or return;
        
    my $writer = Genome::Model::Tools::Consed::Navigation::Writer->new(
        output => $self->acenav,
        title => $reader->title,
    )
        or return;

    while ( my $nav = $reader->next ) {
        $nav->{acefile} = $self->acefile;
        $writer->write_one($nav)
            or return;
    }

    return 1;
}

1;

=pod

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> <ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

