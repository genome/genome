package BAP::DB::TagName;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('tag_name');
__PACKAGE__->columns(All => qw/tag_name/);


1;
