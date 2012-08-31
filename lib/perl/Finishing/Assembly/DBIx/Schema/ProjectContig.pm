package Finishing::Assembly::DBIx::Schema::ProjectContig;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('project_contig');
__PACKAGE__->add_columns
(
    'project_id' => 
    {
      'data_type' => 'int',
      'size' => '11',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
    },
    'contig_id' => 
    {
      'data_type' => 'int',
      'size' => '11',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
    },
);
__PACKAGE__->belongs_to('project', 'Finishing::Assembly::DBIx::Schema::Project', 'project_id');
__PACKAGE__->belongs_to('contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'contig_id');

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/ProjectContig.pm $
#$Id: ProjectContig.pm 30966 2007-12-13 21:28:06Z ebelter $
