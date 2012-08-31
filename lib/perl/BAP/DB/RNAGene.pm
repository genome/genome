package BAP::DB::RNAGene;

use base 'BAP::DB::DBI';


my $start_col      = Class::DBI::Column->new('seq_start'   => { accessor => 'start' });
my $end_col        = Class::DBI::Column->new('seq_end'     => { accessor => 'end' });
my $desc_col       = Class::DBI::Column->new('description' => { accessor => 'desc' });
my $rfam_prod_col  = Class::DBI::Column->new('product'     => { accessor => 'rfam_prod' });


__PACKAGE__->table('rna_gene');
__PACKAGE__->columns(
                     'All' => qw(
                                 gene_id
                                 gene_name
                                 sequence_id
                                 acc
                                 strand
                                 score
                                 source
                                 redundant
                                ),
                     $start_col,
                     $end_col,
                     $desc_col,
		             $rfam_prod_col,
                    );
__PACKAGE__->sequence('gene_id_seq');

__PACKAGE__->has_many('gene_tags' => BAP::DB::GeneTag, 'gene_id');

1;
