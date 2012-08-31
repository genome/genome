package Finishing::Assembly::DBIx::Schema::ConsensusSequence;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('assembly.consensus_sequence');


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
    'contig_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'contig_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'bases' => {
      'data_type' => 'clob',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'bases',
      'is_nullable' => 0,
      'size' => 0,
      bind_attributes => {
          ora_type => 113,
          ora_field => 'bases',
      },
    },
    'qualities' => {
      'data_type' => 'clob',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'qualities',
      'is_nullable' => 0,
      'size' => 0,
      bind_attributes => {
          ora_type => 113,
          ora_field => 'qualities',
      },
    },
);
__PACKAGE__->set_primary_key('id');



__PACKAGE__->belongs_to('contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'contig_id');


1;
