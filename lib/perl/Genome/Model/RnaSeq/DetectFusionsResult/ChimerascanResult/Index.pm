package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult::Index {
    is => 'Genome::Model::RnaSeq::DetectFusionsResult::Index::Base',
    has_param => [
        version => {
            is => 'Text',
            doc => 'the version of chimerascan to use to make the index',
        },
        bowtie_version => {
            is => 'Text',
            doc => 'the version of bowtie to use to make the index',
        },
    ],
    has_input => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'The annotation build from which to derive the gene file',
        },
    ],
    has_calculated => [
        gene_file => {
            is => 'Text',
            is_optional => 1,
            calculate => sub { $_[0] . "/gene_file";},
            calculate_from => [ "temp_staging_directory" ],
        },
    ],

    doc => 'This holds the bowtie indices and modified FASTA required to run chimerascan',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_) or return;

    $self->_prepare_staging_directory;

    $self->prepare_gene_file;
    $self->run_indexer;
    $self->get_sequence_dictionary;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan-index/' . $self->id;
}

#rumour has it future versions of chimerascan will support different formats for this
# TODO: Can we just convert the GTF file that exists for all annotation builds using gtf_to_genepred.py that ships with chimerascan?
sub prepare_gene_file {
    my $self = shift;

    my $annotation_build = $self->annotation_build;

    my $itr = $annotation_build->transcript_iterator;

    my $output_file = Genome::Sys->open_file_for_writing($self->gene_file);

    while (my $tx = $itr->next()){
        next if $tx->coding_region_start eq "NULL";
        next unless $tx->gene_name; #the chimerascan indexer chokes on these

        my @line = ();
        push @line, $tx->transcript_name, $tx->chrom_name;

        #strand needs to be +/- instead of +1/-1
        (my $strand_val = $tx->strand) =~ s/\+1/\+/;
        $strand_val =~ s/-1/-/;

        push @line, $strand_val, $tx->transcript_start, $tx->transcript_stop, $tx->coding_region_start, $tx->coding_region_stop;

        my @start, my @stop = ();
        for (grep {$_->structure_type =~/cds_exon|utr_exon/ } $tx->ordered_sub_structures){
            push @start, $_->structure_start;
            push @stop, $_->structure_stop;
        }

        #yes, these trailing commas are needed. no, I don't like it either.
        push @line, scalar(@start), join(",",@start) . ",", join(",",@stop) . ",", $tx->gene_name;

        $output_file->say(join("\t", @line));

        Genome::Transcript->unload();
        Genome::TranscriptStructure->unload();
    }

    $output_file->close();

    return 1;
}

sub run_indexer {
    my $self = shift;
    my $fasta = $self->reference_build->full_consensus_path('fa');
    my $gene_file = $self->gene_file;
    my $output_dir = $self->temp_staging_directory;

    my $bowtie_dir =  Genome::Model::Tools::Bowtie->base_path($self->bowtie_version);
    my $cmd_path = Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult->_path_for_version($self->version);

    my $cmd = "python $cmd_path/chimerascan_index.py --bowtie-dir=$bowtie_dir $fasta $gene_file $output_dir";

    local $ENV{PYTHONPATH} =  ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "")  .
        Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult->_python_path_for_version($self->version);

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$fasta, $gene_file],

    );

}

sub get_sequence_dictionary {
    my $self = shift;

    my $working_dir;
    if ($self->output_dir && -d $self->output_dir) {
        $working_dir = $self->output_dir;
    } else {
        $working_dir = $self->temp_staging_directory;
    }
    my $fasta_file = $working_dir .'/align_index.fa';
    my $sam_file = $working_dir .'/align_index_tmp.sam';
    my $seqdict_file = $working_dir .'/align_index.sam';

    if (-s $seqdict_file) {
        return $seqdict_file;
    }

    my $species_name = $self->reference_build->species_name;
    my $assembly_name = $self->reference_build->assembly_name;
    unless ($assembly_name) {
        $assembly_name = $self->reference_build->name;
    }
    # gmt picard create-sequence-dictionary
    unless (Genome::Model::Tools::Picard::CreateSequenceDictionary->execute(
        output_file => $sam_file,
        reference_fasta => $fasta_file,
        species => $species_name,
        genome_assembly => $assembly_name,
    )) {
        die('Failed to create sequence dictionary!');
    }

    # gmt picard sort-sam
    unless (Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $sam_file,
        output_file => $seqdict_file,
    )) {
        die('Failed to sort sam file!');
    }
    # TODO:
    die('Validate the order of the chromosome SQ lines with the reference build!');
    unlink($sam_file);
    return $seqdict_file;
}

1;
