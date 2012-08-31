package BAP::DB::GeneTag;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('gene_tag');
__PACKAGE__->columns( All => qw(gene_id tag_id) );

# has a relationships?

__PACKAGE__->has_a('tag_id' => BAP::DB::Tag, );

1;

# $Id$
