package Finishing::Assembly::DBIx::Schema::Library;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('library');


__PACKAGE__->add_columns(
    'prefix' => {
      'data_type' => 'char',
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'prefix',
      'is_nullable' => 0,
      'size' => '5'
    },
    'insert_size' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'insert_size',
      'is_nullable' => 0,
      'size' => '11'
    },
    'insert_stdev' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'insert_stdev',
      'is_nullable' => 0,
      'size' => '11'
    },
);


use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('library');


__PACKAGE__->add_columns(
    'prefix' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'prefix',
      'is_nullable' => 0,
      'size' => '5'
    },
    'insert_size' => {
      'data_type' => 'mediumint',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'insert_size',
      'is_nullable' => 1,
      'size' => '9'
    },
    'insert_stdev' => {
      'data_type' => 'mediumint',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'insert_stdev',
      'is_nullable' => 1,
      'size' => '9'
    },
);
__PACKAGE__->set_primary_key('prefix');

1;
