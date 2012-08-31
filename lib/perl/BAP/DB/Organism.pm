package BAP::DB::Organism;

use base 'BAP::DB::DBI';


__PACKAGE__->table('organism');
__PACKAGE__->columns(
                     'All' => qw(
                                 organism_id
                                 organism_name
                                 ncbi_taxonomy_id
                                 gram_stain
                                 locus
                                )
                     );
__PACKAGE__->sequence('organism_id_seq');

1;
