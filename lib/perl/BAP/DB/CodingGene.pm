package BAP::DB::CodingGene;

use base 'BAP::DB::DBI';
use DBD::Oracle qw(:ora_types);


my $start_col = Class::DBI::Column->new('seq_start' => { accessor => 'start'  });
my $end_col   = Class::DBI::Column->new('seq_end'   => { accessor => 'end'    });

__PACKAGE__->table('mgap.coding_gene');

__PACKAGE__->columns(
                     'All' => qw(
                                 gene_id
                                 gene_name
                                 sequence_id
                                 sequence_string
                                 strand
                                 score
                                 source
                                 internal_stops
                                 missing_start
                                 missing_stop
                                 fragment
                                 wraparound
                                 blastp_evidence
                                 pfam_evidence
                                 phase_0
                                 phase_1
                                 phase_2
                                 phase_3
                                 phase_4
                                 phase_5
                                ), 
                     $start_col,
                     $end_col,
                    );

my $seq_type = { ora_type => ORA_BLOB };

__PACKAGE__->__data_type({});

__PACKAGE__->data_type('sequence_string' => $seq_type);

__PACKAGE__->sequence('gene_id_seq');

__PACKAGE__->has_many('protein' => BAP::DB::Protein, 'gene_id');

__PACKAGE__->has_many('gene_tags' => BAP::DB::GeneTag, 'gene_id');


# not pretty, but sends back the tag names and values with a single call...
sub tag_names
{
    my $self = shift;

    my @tags = $self->gene_tags;
    my @tagreturns;
    foreach my $tag (@tags)
    {
        my $tag_row = $tag->tag_id;
        push(@tagreturns,$tag_row);
    }
    return @tagreturns;
}

1;
