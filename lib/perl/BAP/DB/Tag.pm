package BAP::DB::Tag;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('tag');
__PACKAGE__->columns( All => qw/tag_id tag_name tag_value/ );

__PACKAGE__->has_a('tag_name' => BAP::DB::TagName);
__PACKAGE__->has_a('tag_value' => BAP::DB::TagValue);

1;
