package Finishing::Assembly::DBIx::Schema::Gap;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('gap');
__PACKAGE__->add_columns
(
    id => 
    {
      data_type => 'int',
      is_auto_increment => 1,
      is_nullable => 0,
      size => '11'
    },
    left_contig_id => 
    {
      data_type => 'int',
      is_foreign_key => 1,
      is_nullable => 1,
      size => 11,
    },
    right_contig_id => 
    {
      data_type => 'int',
      is_foreign_key => 1,
      is_nullable => 1,
      size => 11,
    },
    type =>
    {
        data_type => 'varchar',
        is_nullable => 1,
        #default_value => 'contig',
        size => 30,
    },
    linkage =>
    {
        data_type => 'varchar',
        is_nullable => 1,
        #default_value => 'no',
        size => 3,
    },
   'length' => 
    {
      data_type => 'int',
      default_value => 1,
      is_nullable => 1,
      size => 11,
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->belongs_to
(
    'left_contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'left_contig_id'
);
__PACKAGE__->belongs_to
(
    'right_contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'right_contig_id'
);

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/Gap.pm $
#$Id: Gap.pm 30516 2007-11-30 21:38:46Z ebelter $
