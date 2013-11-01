package Genome::Model::Tools::Bacterial::DumpSequences;

use strict;
use warnings;

use Genome;
use Command;
use Carp;
use Bio::Seq;
use Bio::SeqIO;
use BAP::DB::SequenceSet;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'sequence_set_id' => {is => 'Integer',
                              doc => "sequence set id in mgap",  },
        'output_dir' => { is => 'String',
                          doc => "output directory where output files get written", },
        'locus_tag'  => { is => 'String',
                          doc => "locus tag name; used to name output files", },
        'dev'        => { is => 'Boolian',
                          doc => "use development mgap db",
                          is_optional => 1,
                          default => 0,
        },
        'phase'      => { is => 'Integer',
                          doc => "which phase to dump (default: 5)",
                          is_optional => 1,
                          default => 5,
        },
    ],
    );

sub help_brief
{
    "script for dumping sequences from mgap";
}

sub help_detail
{
    return <<EOS
this script dumps protein and cds dna sequence from the MGAP database.
EOS
}

sub execute
{
    my $self = shift;
    # connect
    my $sequence_set = BAP::DB::SequenceSet->retrieve($self->sequence_set_id);
    my @sequences = $sequence_set->sequences;

    my $peptidefile = $self->outputdir ."/". $self->locus_tag .".pep.fa";
    my $cdsfile = $self->outputdir ."/". $self->locus_tag . ".cds.fa";
    my $pepout = Bio::SeqIO->new(-file => ">$peptidefile", -format => 'fasta');
    my $cdsout = Bio::SeqIO->new(-file => ">$cdsfile", -format => 'fasta');
    foreach my $sequence (@sequences)
    {
        # get all the cds sequences
        # next GENE unless phase_5
        my @coding_genes = $sequence->coding_genes;
GENE:        foreach my $cds (@coding_genes) {
            next GENE unless $cds->phase_5;
            my $cds_seq = Bio::Seq->new(-seq => $cds->sequence_string,
                                        -display_id => $cds->gene_name );
            $cdsout->write_seq($cds_seq);
            my @proteins = $cds->protein;
            foreach my $pep (@proteins) 
            {
                my $pep_seq = Bio::Seq->new(-seq => $pep->sequence_string,
                                            -display_id => $cds->protein_name);
                $pepout->write_seq($pep_seq);
            }
            
        }
        # get the protein sequences

        # keep the sequence objects in respective arrays?
        # or do this as we get them?
  
    }
    # grab cds

    # grab proteins

    return 1;
}
