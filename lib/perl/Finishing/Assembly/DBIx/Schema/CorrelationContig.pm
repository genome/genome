package Finishing::Assembly::DBIx::Schema::CorrelationContig;

use base 'DBIx::Class';

use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('correlation_contig');
__PACKAGE__->add_columns(
    'correlation_id' => {
      'data_type' => 'int',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    'contig_id' => {
      'data_type' => 'int',
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
);
__PACKAGE__->belongs_to('correlation', 'Finishing::Assembly::DBIx::Schema::ImprovementCorrelation', 'correlation_id');
__PACKAGE__->belongs_to('contig', 'Finishing::Assembly::DBIx::Schema::Contig', 'contig_id');

1;
