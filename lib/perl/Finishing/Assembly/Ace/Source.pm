package Finishing::Assembly::Ace::Source;

use strict;
use warnings;

use Data::Dumper;

sub new
{
    my ($class, %p) = @_;

    return bless \%p, $class;
}

sub mk_code_accessors
{
    my ($class, @attrs) = @_;
    
    foreach my $attr ( @attrs )
    {
        no strict 'refs';
        *{ $class . '::' . $attr } = sub
        {
            my $self = shift;

            return $self->{'_' . $attr}->(@_);
        };
    }

    return 1;
}

sub DESTROY
{
    my $self = shift;

    $self->{_destroy}->();
};

1;

=pod

=head1 Name

Finishing::Assembly::Ace::Source

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

