package Finishing::Assembly::Commands::CreateOrganism;

use strict;
use warnings;
no warnings 'reserved';

use base 'Finishing::Assembly::Commands::Base';

use Data::Dumper;

my %organism_name :name(organism_name:r)
    :desc('Name of organism to create');

sub execute
{
    my $self = shift;

    my $factory = $self->_factory;

    my $existing_organism = $factory->get_organism( $self->name );
    return $existing_organism if $existing_organism;

    return $factory->create_organism( $self->name );
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::CreateOrganism

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

