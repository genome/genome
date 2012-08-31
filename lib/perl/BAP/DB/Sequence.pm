package BAP::DB::Sequence;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);

__PACKAGE__->table('mgap.dna_sequence');

__PACKAGE__->columns(
                     'All' => qw(
                                 sequence_id
                                 sequence_name
                                 sequence_set_id
                                 sequence_string
                                )
                    );

my $seq_type = { ora_type => ORA_BLOB };

__PACKAGE__->__data_type({});

__PACKAGE__->data_type('sequence_string' => $seq_type);

__PACKAGE__->sequence('sequence_id_seq');

__PACKAGE__->has_a('sequence_set_id' => BAP::DB::SequenceSet);

__PACKAGE__->has_many('coding_genes' => BAP::DB::CodingGene, 'sequence_id');
__PACKAGE__->has_many('trna_genes' => BAP::DB::tRNAGene, 'sequence_id');
__PACKAGE__->has_many('rna_genes' => BAP::DB::RNAGene, 'sequence_id');

1;
