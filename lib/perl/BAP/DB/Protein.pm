package BAP::DB::Protein;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('mgap.protein');

__PACKAGE__->columns(
                     'All' => qw(
                                 protein_id
                                 protein_name
                                 gene_id
                                 sequence_string
                                 internal_stops
                                 cellular_localization 
                                 enzymatic_pathway_id
                                 cog_id
                                )
                    );

my $seq_type = { ora_type => ORA_BLOB };

__PACKAGE__->__data_type({});

__PACKAGE__->data_type('sequence_string' => $seq_type);

__PACKAGE__->sequence('protein_id_seq');

__PACKAGE__->has_a('gene_id' => BAP::DB::CodingGene);

1;
