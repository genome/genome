package Genome::Model::ReferenceSequence::Command::CalculateGenomeTerritory;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Model::ReferenceSequence::Command::CalculateGenomeTerritory {
    is => 'Command::V2',
    has => [
       reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference sequence build to calculate genome territory.',
        },
    ],
    doc => 'Count the number of non-N reference sequence nucleotides.',
};

sub help_detail {
    return <<EOHELP;
This command iterates over all chromosomes of the provided reference sequence and counts the number on non-N nucleotides.  This number represents the total genomic territory represented by this specfic reference sequence build.
EOHELP
}

sub execute {
    my $self = shift;

    my $build = $self->reference_sequence_build;
    my $metric_name = 'GENOME_TERRITORY';
    my $metric_value = $build->get_metric($metric_name);
    if ( defined($metric_value) ) {
        $self->error_message($metric_name .' already exists for build!');
        return;
    }
    my $fasta_file = $build->full_consensus_path('fa');
    my $total_bases = `grep -v '>' $fasta_file | tr -cd [:alpha:] | wc -c`;
    my $n_count = `grep -v '>' $fasta_file | tr -cd [Nn] | wc -c`;
    $metric_value = $total_bases - $n_count;
    $build->set_metric($metric_name,$metric_value);

    return 1;
}


1;
