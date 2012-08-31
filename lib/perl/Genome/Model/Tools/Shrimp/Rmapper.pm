package Genome::Model::Tools::Shrimp::Rmapper;

class Genome::Model::Tools::Shrimp::Rmapper {
    is => 'Genome::Model::Tools::Shrimp',
    has_input => [
        fasta_file => {
            is => 'Text',
            doc => 'The reads in a fasta format file',
        },
        reference_directory => {
            is => 'Text',
            doc => 'The reference genome directory where chromosome fasta files exist',
        },
        output_directory => {
            is => 'Text',
            doc => 'The base directory where all output is written',
        },
        aligner_params => {
            is => 'Text',
            is_optional => 1,
            default_value => '',
            doc => 'A string of aligner params and options.  See SHRiMP documentation for details.'
        },
    ],
    has_output  => [
        aligner_output_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The file where shrimp aligner messages are output',
        },
        alignment_file => {
            is => 'Text',
            is_optional => 1,
            doc => 'The file where shrimp alignments are output',
        },
    ],
};

sub execute {
    my $self = shift;

    my $basename = File::Basename::basename($self->fasta_file);
    my @ref_fasta_files = grep { $_ !~ /all_sequences/ } glob($self->reference_directory .'/*.fasta');
    $self->alignment_file($self->output_directory .'/'. $basename .'.shrimp');
    $self->aligner_output_file($self->output_directory .'/'. $basename .'.aligner_output');
    my $params = $self->aligner_params;
    my $cmd = sprintf('%s %s %s %s 1> %s 2>> %s',
                      $self->shrimp_path .'/rmapper-'.$self->read_space,
                      $params,
                      $self->fasta_file,
                      join(' ',@ref_fasta_files),
                      $self->alignment_file,
                      $self->aligner_output_file);

    Genome::Sys->shellcmd(
        cmd          => $cmd,
        input_files  => [ $self->fasta_file, @ref_fasta_files ],
        output_files => [ $self->alignment_file, $self->aligner_output_file ],
        skip_if_output_is_present => 0,
    );
    return 1;
}

1;
