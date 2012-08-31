package BAP::DB::tRNAGene;

use base 'BAP::DB::DBI';


my $start_col = Class::DBI::Column->new('seq_start' => { accessor => 'start' });
my $end_col   = Class::DBI::Column->new('seq_end'   => { accessor => 'end'   });

__PACKAGE__->table('trna_gene');
__PACKAGE__->columns(
                     'All' => qw(
                                 gene_id
                                 gene_name
                                 sequence_id
                                 codon
                                 aa
                                 strand
                                 score
                                 source
                                ),
                     $start_col, 
                     $end_col,
                    );
__PACKAGE__->sequence('gene_id_seq');

__PACKAGE__->has_many('gene_tags' => BAP::DB::GeneTag , 'gene_id');

1;
