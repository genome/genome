package Genome::Utility::Window::Variation;

use strict;
use warnings;

class Genome::Utility::Window::Variation{
    is => 'Genome::Utility::Window'
};

use Data::Dumper;

sub object_start_method
{
    return 'start';
}

sub object_stop_method
{
    return 'stop';
}

sub variations
{
    my $self = shift;

    return @{ $self->_objects };
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

