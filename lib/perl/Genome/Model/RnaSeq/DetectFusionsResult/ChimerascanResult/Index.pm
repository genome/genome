package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanBase::Index',
};

sub run_indexer {
    my $self = shift;
    my $fasta = $self->reference_build->full_consensus_path('fa');
    my $gene_file = $self->gene_file;
    my $output_dir = $self->temp_staging_directory;

    my $bowtie_dir =  Genome::Model::Tools::Bowtie->base_path($self->bowtie_version);
    my $cmd_path = Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult->_path_for_version($self->version);

    my $cmd = "python $cmd_path/chimerascan_index.py --bowtie-dir='$bowtie_dir' '$fasta $gene_file' '$output_dir'";

    local $ENV{PYTHONPATH} =  ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "")  .
        Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult->_python_path_for_version($self->version);

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta, $gene_file],
    );
}

1;
