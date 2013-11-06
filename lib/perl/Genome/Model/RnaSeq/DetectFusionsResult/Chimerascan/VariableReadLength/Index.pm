package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Index;

use strict;
use warnings;

use Genome;
use File::Spec qw();

class Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::IndexBase',
};

sub run_indexer {
    my $self = shift;
    my $fasta = $self->reference_build->full_consensus_path('fa');
    my $gene_file = $self->gene_file;
    my $output_dir = $self->temp_staging_directory;

    my $bowtie_dir =  Genome::Model::Tools::Bowtie->base_path($self->bowtie_version);
    my $executable = Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result->get_executable_path($self->version);

    $self->_run($executable, $bowtie_dir, $fasta, $gene_file, $output_dir);
}

sub _run {
    my ($self, $executable, $bowtie_dir, $fasta, $gene_file, $output_dir) = @_;

    my $cmd = "$executable chimerascan_index.py --bowtie-dir='$bowtie_dir' '$fasta' '$gene_file' '$output_dir'";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$executable, $fasta, $gene_file],
        output_directories => [$output_dir],
        redirect_stderr => File::Spec->join($output_dir, 'chimerascan_index.stderr'),
    );
}

sub prepare_gene_file {
    my $self = shift;

    my $annotation_build = $self->annotation_build;
    my $reference_sequence = $annotation_build->reference_sequence;
    die("Couldn't find reference_sequence on annotation_build (" . $annotation_build->id . ").") unless $reference_sequence;

    my $gtf_file = $annotation_build->annotation_file('gtf', $reference_sequence->id);
    my $gene_file = $self->gene_file;
    $self->_convert_gtf_to_genepred($gtf_file, $gene_file);
    return 1;
}

sub _convert_gtf_to_genepred {
    my ($self, $gtf_file, $gene_file) = @_;

    my $executable = Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result->get_executable_path($self->version);
    my $cmd = "$executable gtf_to_genepred.py '$gtf_file' '$gene_file'";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$gtf_file],
        output_files => [$gene_file],
        redirect_stdout => '/dev/null',
        redirect_stderr => '/dev/null',
    );
    return 1;
}

1;
