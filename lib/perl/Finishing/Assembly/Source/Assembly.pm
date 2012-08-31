package Finishing::Assembly::Source::Assembly;

use strict;
use warnings;

use base 'Finfo::Accessor';

use Data::Dumper;

__PACKAGE__->mk_accessors(qw/ contigs /);

sub get_contig
{
    my ($self, $name) = @_;

    return unless $self->contigs;

    foreach my $ctg ( @{ $self->contigs } )
    {
        return $ctg if $ctg->name eq $name;
    }

    return;
}

sub get_assembled_read
{
    my ($self, $name) = @_;

    return unless $self->contigs;

    foreach my $ctg ( @{ $self->contigs } )
    {
        my $ri = $ctg->assembled_reads;
        # TODO
    }

    return;
}

#- TAGS -#

sub tags{
    my ($self, $tags) = @_;
    $self->{tags} = [] unless $self->{tags};
    $self->{tags} = $tags if $tags;
    return $self->{tags}; 
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

