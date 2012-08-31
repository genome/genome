package Finishing::Assembly::Phd::AssembledRead;

use strict;
use warnings;

use base 'Finfo::Accessor';

use Data::Dumper;

__PACKAGE__->mk_accessors
(qw/
    name
    base_string qualities
    time chromat_file chem dye phd_file
    /);

sub tags
{
    my ($self, $new_tags) = @_;

    $self->{tags} = $new_tags if $new_tags;

    return $self->{tags} || [];
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::AssembledRead

=head1 Synopsis

=head1 Usage

=head1 Methods

=head1 See Also

=head1 Disclaimer

Copyright (C) 2007 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

