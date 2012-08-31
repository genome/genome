package Genome::Model::Build::DeNovoAssembly::Abyss;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';

class Genome::Model::Build::DeNovoAssembly::Abyss {
    is => 'Genome::Model::Build::DeNovoAssembly',
};

sub fastq_input_files {
    my $self = shift;
    my $dir = $self->data_directory;
    return ("$dir/fwd.fq", "$dir/rev.fq");
}
    

sub existing_assembler_input_files {
    my $self = shift;
    return grep { -s $_ } $self->fastq_input_files;
}

sub read_processor_output_files_for_instrument_data {
    return $_[0]->fastq_input_files;
}

sub contigs_fasta_file {
    return $_[0]->data_directory.'/all-contigs.fa';
}

# assemble
sub assembler_params {
    my $self = shift;

    my %params = $self->processing_profile->assembler_params_as_hash;
    $params{version} = $self->processing_profile->assembler_version;

    ($params{fastq_a}, $params{fastq_b}) = $self->fastq_input_files;

    my $output_dir = $self->data_directory;
    $params{output_directory} = $output_dir;

    return %params;
}

1;
