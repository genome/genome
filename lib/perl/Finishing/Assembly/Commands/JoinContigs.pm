package Finishing::Assembly::Commands::JoinContigs;

use strict;
use warnings;

use base 'Finishing::Assembly::Commands::AnalyzeJoin';

use Data::Dumper;

sub execute
{
    my $self = shift;

    my $contig_tools = Finishing::Assembly::ContigTools->new;
    my $new_contig = $contig_tools->merge
    (
        $self->_left_contig,
        $self->_right_contig,
        $self->_phd, #TODO
        %{ $self->_contig_tools_params },
    );

    #TODO set new contig
    
    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::JoinContigs

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
