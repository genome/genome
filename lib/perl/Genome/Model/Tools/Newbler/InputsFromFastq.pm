package Genome::Model::Tools::Newbler::InputsFromFastq;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Newbler::InputsFromFastq {
    is => 'Genome::Model::Tools::Newbler',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Path to assembly',
        },
    ],
};

sub help_brief {
    'Tool to create input fasta and qual files from input fastq';
}

sub help_detail {
    return <<"EOS"
gmt newbler inputs-from-fastq --assembly-directory /gscmnt/999/Newbler_e_coli_assembly
EOS
}

sub execute {
    my $self = shift;

    unless ( -d $self->consed_edit_dir ) {
        $self->create_consed_dir;
    }

    my @fastq_files = $self->input_fastq_files;
    for my $fastq_file ( @fastq_files ) {
        my $root_name = File::Basename::basename( $fastq_file );
        $root_name =~ s/\.fastq$//;
        my $fasta_out = $self->assembly_directory."/consed/edit_dir/$root_name".'.fasta';
        my $qual_out = $self->assembly_directory."/consed/edit_dir/$root_name".'.fasta.qual';
        my $fq = Genome::Model::Tools::Sx->create(
            input => [ $fastq_file ],
            output => [ $fasta_out.':qual_file='.$qual_out ],
        );
        if ( not $fq ) {
            $self->error_message( "Failed to create fastq to fasta converter" );
            return;
        }
        if ( not $fq->execute ) {
            $self->error_message( "Failed to execute fastq to fasta converter" );
            return;
        }
    }

    return 1;
}

1;
