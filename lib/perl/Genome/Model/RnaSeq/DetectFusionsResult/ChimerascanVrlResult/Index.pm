package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanVrlResult::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanVrlResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanBase::Index',
};

sub run_indexer {
    my $self = shift;
    my $fasta = $self->reference_build->full_consensus_path('fa');
    my $gene_file = $self->gene_file;
    my $output_dir = $self->temp_staging_directory;

    my $bowtie_dir =  Genome::Model::Tools::Bowtie->base_path($self->bowtie_version);
    my $executable = Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanVrlResult->get_executable_path($self->version);

    my $cmd = "$executable chimerascan_index.py --bowtie-dir=$bowtie_dir $fasta $gene_file $output_dir";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta, $gene_file],
    );
}

1;
