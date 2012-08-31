package Finishing::Assembly::DBIx::Schema::ConsensusTag;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core InflateColumn::DateTime/);
__PACKAGE__->table('consensus_tag');


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
    'no_trans' => {
      'data_type' => 'boolean',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'no_trans',
      'is_nullable' => 1,
      'size' => 0
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
      'size' => '40'
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



__PACKAGE__->belongs_to('contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'contig_id');

#TODO add oligo and autoexp accessor methods (implement following interface)

    my %oligo_methods = (
        oligo_name  => undef,
        oligo_seq  => undef,
        oligo_temp  => undef,
        oligo_templates  => undef,
        complemented => undef,
        orientation => undef,
    );

    my %auto_exp_methods = (

        orientation => undef,
        num1 => undef,
        num2 => undef,
        num3 => undef,
        chem => undef,
        primer_type => undef,
        purpose  => undef,
        fix_cons_errors => undef,
        original_cons_errors => undef,
        original_single_subclone_bases => undef, 
        primer => undef,
        temp => undef,
        id => undef,
        exp_id_and_template => undef,
        oligo_name => undef,
        oligo_seq => undef,
        oligo_temp => undef,
    );

1;
