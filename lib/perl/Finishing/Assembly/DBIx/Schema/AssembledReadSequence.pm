package Finishing::Assembly::DBIx::Schema::AssembledReadSequence;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('assembled_read_sequence');


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
    'assembled_read_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'assembled_read_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'bases' => {
      'data_type' => 'blob',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'bases',
      'is_nullable' => 0,
      bind_attributes => {
          ora_type => 113,
          ora_field => 'bases',
      },
      'size' => 0
    },
    'qualities' => {
      'data_type' => 'blob',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'qualities',
      'is_nullable' => 0,
#      'size' => '4294967295'
      bind_attributes => {
          ora_type => 113,
          ora_field => 'qualities',
      },
      'size' => 0
    },
);
__PACKAGE__->set_primary_key('id');



__PACKAGE__->belongs_to('assembled_read', 'Finishing::Assembly::DBIx::Schema::AssembledRead', 'assembled_read_id');


1;
