package Finishing::Assembly::Ace::AssembledRead;

use strict;
use warnings;

use base 'Finishing::Assembly::Ace::Source';

use Data::Dumper;

__PACKAGE__->mk_code_accessors
(qw/ 
    name rename 
    length position complemented
    base_string 
    tags info_count
    qual_clip_start qual_clip_stop align_clip_start align_clip_stop
    time chromat_file chem dye phd_file
    /);

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

