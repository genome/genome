package Finishing::Assembly::DBIx::Schema::ReplacedContig;

use strict;
use warnings;

use base 'DBIx::Class';

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('replaced_contig');
#    id integer not null primary key,
#    replacement_id integer not null
__PACKAGE__->add_columns(
    'id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    'replacement_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
);
__PACKAGE__->set_primary_key('id');
__PACKAGE__->has_one
(
    'replaced_contig_event',
    'Finishing::Assembly::DBIx::Schema::ReplacedContigEvent',
    'old_contig_id'
);
__PACKAGE__->might_have
(
    'replacing_contig_event',
    'Finishing::Assembly::DBIx::Schema::ReplacedContigEvent',
    'new_contig_id'
);

1;
