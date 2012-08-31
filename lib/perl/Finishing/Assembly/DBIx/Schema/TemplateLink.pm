package Finishing::Assembly::DBIx::Schema::TemplateLink;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('template_link');
__PACKAGE__->add_columns(
    'template_name' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'template_name',
      'is_nullable' => 0,
      'size' => '50'
    },
    'left_contig_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'left_contig_id',
      'is_nullable' => 1,
      'size' => '11'
    },
    'right_contig_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'right_contig_id',
      'is_nullable' => 1,
      'size' => '11'
    },
    'left_read_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'left_read_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'left_read_position' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'left_read_position',
      'is_nullable' => 0,
      'size' => '11'
    },
    'right_read_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'right_read_id',
      'is_nullable' => 0,
      'size' => '11'
    },
    'right_read_position' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'right_read_position',
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
    
);



__PACKAGE__->belongs_to('left_contig','Finishing::Assembly::DBIx::Schema::Contig',
    'left_contig_id',
);

__PACKAGE__->belongs_to('right_contig','Finishing::Assembly::DBIx::Schema::Contig',
    'right_contig_id'
);
__PACKAGE__->belongs_to('left_read','Finishing::Assembly::DBIx::Schema::AssembledRead',
    'left_read_id'
);
__PACKAGE__->belongs_to('right_read','Finishing::Assembly::DBIx::Schema::AssembledRead',
    'right_read_id'
);
__PACKAGE__->belongs_to('assembly','Finishing::Assembly::DBIx::Schema::Assembly',
    'assembly_id'
);

sub left_scaffold {
    my $self = shift;
    return $self->left_contig->scaffold;
}

sub right_scaffold {
    my $self = shift;
    return $self->right_contig->scaffold;
}

1;
