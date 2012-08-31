package Finishing::Assembly::DBIx::Schema::ChromosomeScaffoldPosition;
# TODO remove this package and table

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('chromosome_scaffold_position');
__PACKAGE__->add_columns
(
    chromosome_id =>
    {
      data_type => 'int',
      is_auto_increment => 0,
      default_value => undef,
      is_foreign_key => 1,
      is_nullable => 0,
      size => 11,
    },
    scaffold_id =>
    {
      data_type => 'int',
      is_auto_increment => 0,
      default_value => undef,
      is_foreign_key => 1,
      is_nullable => 0,
      size => 11,
    },
    start_position => 
    {
      data_type => 'int',
      is_auto_increment => 0,
      default_value => undef,
      is_nullable => 1,
      size => 11,
    },
    stop_position => 
    {
      data_type => 'int',
      is_auto_increment => 0,
      default_value => undef,
      is_nullable => 1,
      size => 11,
    },
);

__PACKAGE__->set_primary_key('chromosome_id', 'scaffold_id');
__PACKAGE__->add_unique_constraint([qw/ chromosome_id scaffold_id /]);
__PACKAGE__->belongs_to('scaffold', 'Finishing::Assembly::DBIx::Schema::Scaffold', 'scaffold_id');
__PACKAGE__->belongs_to('chromosome', 'Finishing::Assembly::DBIx::Schema::Chromosome', 'chromosome_id');

sub orientation
{
    my $self = shift;

    return '+' unless defined $self->start_position;

    return ( $self->start_position <= $self->stop_position ) ? '+' : '-';
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/ChromosomeScaffoldPosition.pm $
#$Id: ChromosomeScaffoldPosition.pm 30516 2007-11-30 21:38:46Z ebelter $
