package Finishing::Assembly::DBIx::Schema::AssembledRead;

use base 'DBIx::Class';

use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('assembled_read');


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
    'name' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'name',
      'is_nullable' => 0,
      'size' => '50'
    },
    'template_name' => {
      'data_type' => 'varchar',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'template_name',
      'is_nullable' => 0,
      'size' => '50'
    },
    'direction' => {
      'data_type' => 'char',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'direction',
      'is_nullable' => 0,
      'size' => '1'
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
    'contig_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 1,
      'name' => 'contig_id',
      'is_nullable' => 1,
      'size' => '11'
    },
    'length' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'length',
      'is_nullable' => 1,
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
    'complemented' => {
      'data_type' => 'int',
      'is_auto_increment' => 0,
      'default_value' => undef,
      'is_foreign_key' => 0,
      'name' => 'complemented',
      'is_nullable' => 1,
      'size' => '11'
    },
);
__PACKAGE__->set_primary_key('id');


__PACKAGE__->belongs_to('assembly', 'Finishing::Assembly::DBIx::Schema::Assembly','assembly_id');
__PACKAGE__->belongs_to('contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'contig_id',
    {
        proxy => [qw/scaffold/],
    }
);

__PACKAGE__->belongs_to(
    'mate', 
    'Finishing::Assembly::DBIx::Schema::AssembledRead',
    {
        'foreign.template_name' => 'self.template_name',
        'foreign.assembly_id'   => 'self.assembly_id',
    },
    {
        where => {
            'me.id' => \"<> mate.id",
        },
        join => 'mate',
    },
);

sub scaffold_position {
    my $self = shift;
    my $contig = $self->contig;
    return unless $contig;
    return $contig->scaffold_position + $self->start_position;
}

sub library {
    my $self = shift;
    my $name = $self->name;
    my $prefix = substr($name, 0, 5);
    my $schema = $self->result_source->schema;
    my $library = $schema->resultset('Library')->find(
        {
            prefix => $prefix,
        },
    );
    return $library;
}

sub insert_size {
    my $self = shift;
    my $library = $self->library;
    if ($library) {
        return $library->insert_size;
    }
    return undef;
}

sub insert_stdev {
    my $self = shift;
    my $library = $self->library;
    if ($library) {
        return $library->insert_stdev;
    }
    return undef;
}


__PACKAGE__->has_many('tags', 'Finishing::Assembly::DBIx::Schema::ReadTag', 'read_id');
__PACKAGE__->has_one('sequence', 'Finishing::Assembly::DBIx::Schema::AssembledReadSequence', 'assembled_read_id');

=cut
__PACKAGE__->add_relationship('mate', 'Finishing::Assembly::DBIx::Schema::AssembledRead', 
    {
        'foreign.template_name' => 'self.template_name',
        'foreign.assembly_id' => 'self.assembly_id',
    },
);
=cut

1;
