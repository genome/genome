package Finishing::Assembly::DBIx::Schema::AssemblyTag;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core InflateColumn::DateTime/);
__PACKAGE__->table('assembly_tag');


__PACKAGE__->add_columns(
    'id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'assembly_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'assembly_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'tag_type' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'tag_type',
      'is_nullable' => 0,
      'size' => '50'
    },
    'source' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'source',
      'is_nullable' => 1,
      'size' => '20'
    },
    'creation_date' => {
      'data_type' => 'datetime',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'creation_date',
      'is_nullable' => 0,
      'size' => 0
    },
    'text' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'data',
      'is_nullable' => 1,
      'size' => '500'
    },
);
__PACKAGE__->set_primary_key('id');


__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly', 'assembly_id');


1;
