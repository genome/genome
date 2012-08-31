package Finishing::Assembly::DBIx::Schema::ReadTag;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core InflateColumn::DateTime/);
__PACKAGE__->table('read_tag');


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
    'read_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'read_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'start_position' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'start_position',
      'is_nullable' => 1,
      'size' => '11'
    },
    'stop_position' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'stop_position',
      'is_nullable' => 1,
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
    'tag_comment' => {
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


__PACKAGE__->belongs_to('read', 'Finishing::Assembly::DBIx::Schema::AssembledRead', 'read_id');


1;
