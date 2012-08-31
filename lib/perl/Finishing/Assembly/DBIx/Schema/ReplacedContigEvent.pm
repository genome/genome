package Finishing::Assembly::DBIx::Schema::ReplacedContigEvent;
use base 'DBIx::Class';
use strict;
use warnings;

__PACKAGE__->load_components(qw/ Core/);
__PACKAGE__->table('replaced_contig_event');
__PACKAGE__->add_columns
(
    old_contig_id => 
    {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
    new_contig_id => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_foreign_key' => 0,
      'is_nullable' => 0,
      'size' => '11'
    },
    'event_id' => {
      'data_type' => 'int',
      'is_auto_increment' => 1,
      'is_foreign_key' => 1,
      'is_nullable' => 0,
      'size' => '11'
    },
);
__PACKAGE__->belongs_to('event', 'Finishing::Assembly::DBIx::Schema::Event', 'event_id');
__PACKAGE__->belongs_to('old_contig', 'Finishing::Assembly::DBIx::Schema::ReplacedContig', 'old_contig_id');
__PACKAGE__->might_have
(
    'new_contig', 
    'Finishing::Assembly::DBIx::Schema::Contig',
    {'foreign.id' => 'self.new_contig_id'}
);
__PACKAGE__->might_have
(
    'new_contig_replaced', 
    'Finishing::Assembly::DBIx::Schema::ReplacedContig', 
    {'foreign.id' =>'self.new_contig_id'}
);

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/branches/adukes/AssemblyRefactor/Schema/ReplacedContigEvent.pm $
#$Id: ReplacedContigEvent.pm 31279 2007-12-24 17:16:03Z ebelter $#
