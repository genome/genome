package Finishing::Assembly::Ace::Assembly;

use strict;
use warnings;

use base 'Finishing::Assembly::Ace::Source';

use Data::Dumper;

__PACKAGE__->mk_code_accessors(qw/
    contig_count assembled_read_count
    contig_names contigs get_contig
    assembled_reads get_assembled_read
    tags
    /);

sub read_count
{
    return assembled_read_count(@_);
}

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Assembly

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

