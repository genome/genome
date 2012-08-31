package Finishing::Assembly::DBIx::Schema::ChromosomeFirstScaffold;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('chromosome_first_scaffold');
__PACKAGE__->add_columns
(
    chromosome_id =>
    {
      data_type => 'int',
      size => 11,
      default_value => undef,
      is_foreign_key => 1,
      is_nullable => 0,
    },
    assembly_id =>
    {
      data_type => 'int',
      size => 11,
      default_value => undef,
      is_foreign_key => 1,
      is_nullable => 0,
    },
    scaffold_id =>
    {
      data_type => 'int',
      size => 11,
      default_value => undef,
      is_foreign_key => 1,
      is_nullable => 0,
    },
);
__PACKAGE__->set_primary_key('chromosome_id', 'assembly_id');
__PACKAGE__->belongs_to('chromosome', 'Finishing::Assembly::DBIx::Schema::Chromosome', 'chromosome_id');
__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly', 'assembly_id');
__PACKAGE__->belongs_to('scaffold', 'Finishing::Assembly::DBIx::Schema::Scaffold', 'scaffold_id');

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/ChromosomeFirstScaffold.pm $
#$Id: ChromosomeFirstScaffold.pm 30649 2007-12-06 19:45:30Z ebelter $
