package Finishing::Assembly::Commands::Base;

use strict;
use warnings;
no warnings 'reserved';

use Finfo::Std;

use Data::Dumper;
use Finishing::Assembly::Factory;

my %db :name(db:r)
    :isa(['in_list', Finishing::Assembly::Factory->available_dbs ])
    :desc('Database to connect to');
my %db_file :name(db_file:o)
    :isa('file_rw')
    :desc('File to use when connecting to database ');
my %factory :name(_factory:p);


sub START : CUMULATIVE
{
    my $self = shift;
    
    return $self->_factory( Finishing::Assembly::Factory->connect($self->db, $self->db_file) );
}

sub DEMOLISH
{
    my $self = shift;

    #TODO
    #$self->_factory->commit;
    #$self->_factory->disconnect;
    
    return 1;
}

1;

=pod

=head1 Name

Finishing::Assembly::Commands::Base

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

