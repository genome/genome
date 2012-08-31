package Finishing::Assembly::DBIx::Schema::Project;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core /);
__PACKAGE__->table('project');
__PACKAGE__->add_columns
(
    id => 
    {
        data_type => 'int',
        size => 11,
        is_auto_increment => 1,
        is_nullable => 0,
    },
    name => 
    {
        data_type => 'varchar',
        size => 20,
        is_nullable => 0,
    },
    directory => 
    {
        data_type => 'varchar',
        size => 128,
        is_nullable => 0,
    },
    organism_id => 
    {
        data_type => 'int',
        size => 11,
        is_nullable => 0,
    },
    project_type => 
    {
        data_type => 'varchar',
        size => 20,
        is_nullable => 1,
        default_value => 'finishing',
    },
    priority =>
    {
        data_type => 'int',
        size => 2,
        is_nullable => 1,
        default_value => 0,
    },
    owner => 
    {
        data_type => 'varchar',
        size => 10,
        is_nullable => 1,
    },
    comments => 
    {
        data_type => 'varchar',
        size => 100,
        is_nullable => 1,
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->belongs_to('organism', 'Finishing::Assembly::DBIx::Schema::Organism', 'organism_id');
__PACKAGE__->has_many('project_contigs', 'Finishing::Assembly::DBIx::Schema::ProjectContig', 'project_id');
__PACKAGE__->many_to_many('contigs', 'project_contigs', 'contig');

1;

#$HeadURL$
#$Id$
