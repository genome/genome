package Genome::Model::Tools::Bacterial::ExportTranscripts;


use strict;
use warnings;
use Genome;

use BAP::DB::SequenceSet;
use BAP::DB::Sequence;

use Bio::Seq;
use Bio::SeqIO;



UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
    has => [
        sequence_set_id => { is => 'Integer',
                             doc => "sequence set id of genes to dump out",
                           },
    ],
    has_optional => [
        phase => {
            is => 'Integer',
            doc => "specify which phase of gene merging to dump from",
            default => 5,
                 },
        dev => {
            is => 'Boolean',
            doc => "use development database",
            default => 0,
               },

    ],
);


sub help_brief {
"Used to dump transcripts for genes from BAP/MGAP database for a certain phase."
}

sub help_detail {
return <<EOS
This script is for dumping the gene transcript sequence from the MGAP database.
EOS

}

sub help_synopsis {
return <<EOS
gmt bacterial export-transcripts --sequence-set-id <ssid> [--dev]
EOS
}


sub execute {
    my $self = shift;
    my $phase = $self->phase;
    my $sequence_set_id = $self->sequence_set_id;
    
    $phase = "phase_$phase";

    if ($self->dev) { $BAP::DB::DBI::db_env = 'dev'; }

    my $fasta_out = Bio::SeqIO->new(-format => 'Fasta', -fh => \*STDOUT);

    my $sequence_set = BAP::DB::SequenceSet->retrieve($sequence_set_id);

    my @sequences = $sequence_set->sequences();

    foreach my $sequence (@sequences) {

        my @coding_genes = $sequence->coding_genes($phase => 1);

        foreach my $coding_gene (@coding_genes) {
            my $seq_obj = Bio::Seq->new(-seq => $coding_gene->sequence_string,
                                        -display_id => $coding_gene->gene_name);

            $fasta_out->write_seq($seq_obj);
        }
    }

    return 1;
}


1;
