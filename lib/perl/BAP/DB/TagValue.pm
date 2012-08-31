package BAP::DB::TagValue;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('tag_value');
__PACKAGE__->columns(All => qw/tag_value/);


1;
