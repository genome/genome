package Genome::Model::Build::DeNovoAssembly::Velvet;

use strict;
use warnings;

use Genome;

require Carp;
use Data::Dumper 'Dumper';

class Genome::Model::Build::DeNovoAssembly::Velvet {
    is => 'Genome::Model::Build::DeNovoAssembly',
};

#< Files >#
sub existing_assembler_input_files {
    my $self = shift;

    my $collated_fastq_file = $self->collated_fastq_file;
    return $collated_fastq_file if -s $collated_fastq_file;

    return;
}

sub collated_fastq_file {
    return $_[0]->data_directory.'/collated.fastq';
}

sub read_processor_output_files_for_instrument_data {
    return $_[0]->collated_fastq_file;
}

sub assembly_afg_file {
    return $_[0]->data_directory.'/velvet_asm.afg';
}

sub contigs_fasta_file {
    return $_[0]->data_directory.'/contigs.fa';
}

sub sequences_file {
    return $_[0]->data_directory.'/Sequences';
}

sub velvet_fastq_file {
    return $_[0]->data_directory.'/velvet.fastq';
}

sub velvet_ace_file {
    return $_[0]->data_directory.'/edit_dir/velvet_asm.ace';
}

sub ace_file {
    return $_[0]->edit_dir.'/velvet_asm.ace';
}

sub _additional_metrics {
    my ($self, $metrics) = @_;

    my $kmer_used = $self->assembler_kmer_used;
    return if not defined $kmer_used;

    $metrics->{assembler_kmer_used} = $kmer_used;

    return 1;
}

sub assembler_params {
    my $self = shift;

    my %params = $self->processing_profile->assembler_params_as_hash;
    $params{version} = $self->processing_profile->assembler_version;

    # die if these params are specified
    for my $calculated_param (qw/genome_len ins_length/) {
        if ( exists $params{$calculated_param} ) {
            Carp::confess (
                $self->error_message("Can't specify $calculated_param as assembler_param .. it will be derived from input data")
            );
        }
    }

    # params that need to be cleaned up
    if ( defined $params{hash_sizes} ) { 
        $params{hash_sizes} = [ split(/\s+/, $params{hash_sizes}) ],
    }

    #params that need to be derived
    my $collated_fastq_file = $self->collated_fastq_file;
    $params{file} = $collated_fastq_file;

    my $genome_len = $self->genome_size;
    $params{genome_len} = $genome_len;

    my $ins_length = $self->calculate_average_insert_size;
    $params{ins_length} = $ins_length if defined $ins_length;

    my $output_dir = $self->data_directory;
    $params{output_dir} = $output_dir;

    return %params;
}

sub after_assemble {
    my $self = shift;

    # contigs fasta files
    my @contigs_fastas_to_remove = glob($self->data_directory.'/*contigs.fa');
    unless ( @contigs_fastas_to_remove ) { # error here??
        $self->error_message("No contigs fasta files produced from running one button velvet.");
        return;
    }

    my $final_contigs_fasta = $self->contigs_fasta_file;
    for my $contigs_fasta_to_remove ( @contigs_fastas_to_remove ) {
        next if $contigs_fasta_to_remove eq $final_contigs_fasta;
        unless ( unlink $contigs_fasta_to_remove ) {
            $self->error_message(
                "Can't remove unnecessary contigs fasta ($contigs_fasta_to_remove): $!"
            );
            return;
        }
    }

    for my $glob (qw/ logfile timing /) {
        for my $file ( glob($self->data_directory.'/*-'.$glob) ) {
            unless ( unlink $file ) {
                $self->error_message("Can't remove unnecessary file ($glob => $file): $!");
                return;
            }
        }
    }

    return 1;
}
#</ASSEMBLE>#

#< Kmer Used in Assembly >#
sub assembler_kmer_used {
    my $self = shift;

    my $velvet_log = $self->data_directory.'/Log';
    my $fh = eval{ Genome::Sys->open_file_for_reading($velvet_log); };
    if ( not $fh ) {
        $self->error_message("Cannot open velvet log file ($velvet_log) to get assembler kmer used.");
        return;
    }

    $fh->getline;
    my $line = $fh->getline;
    return if not $line;

    $line =~ s/^\s+//;
    my @tokens = split(/\s+/, $line);

    return $tokens[2];
}

#< DIFF >#
sub files_ignored_by_diff {
    my $self = shift;
    my @ignored = $self->SUPER::files_ignored_by_diff;
    push @ignored, qw/ Log scaffolds.stor /;
    return @ignored;
}

sub dirs_ignored_by_diff {
    return qw/ logs reports edit_dir /;
}

#TODO - it should test stats.txt and contigs.fa files but this will error since test method
#thinks contigs.fa and supercontigs.fasta are multiple versions of the same file because of the
#way it's grepping for the files .. stats.txt file exists in two places and test method does not
#like that
sub regex_files_for_diff { 
    return qw/ Log Sequences build.xml collated.fastq velvet_asm.afg /;
}
#</ DIFF >#

sub resolve_assemble_lsf_resource {
    return "-n 4 -R 'span[hosts=1] select[type==LINUX64 && mem>30000] rusage[mem=30000]' -M 30000000";
}

1;
