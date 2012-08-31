package Finishing::Assembly::DBIx::Schema::Event;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('event');

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
    'type' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'type',
      'is_nullable' => 0,
      'size' => '50'
    },
    'created_by' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'created_by',
      'is_nullable' => 0,
      'size' => '10'
    },
    'created' => {
      'data_type' => 'timestamp',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'created',
      'is_nullable' => 1,
      'size' => 0
    },
);
__PACKAGE__->set_primary_key('id');

1;
