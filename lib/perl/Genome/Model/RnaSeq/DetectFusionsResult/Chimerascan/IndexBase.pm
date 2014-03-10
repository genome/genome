package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::IndexBase;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::IndexBase {
    is => 'Genome::SoftwareResult::Stageable',
    has_param => [
        version => {
            is => 'Text',
            doc => 'the version of chimerascan to use to make the index',
        },
        bowtie_version => {
            is => 'Text',
            doc => 'the version of bowtie to use to make the index',
        },
        picard_version => {
            is => 'Text',
            doc => 'the version of picard used to manipulate BAM files',
        },
    ],
    has_input => [
        reference_build => {
            is => "Genome::Model::Build::ReferenceSequence",
            doc => 'object representing the reference sequence version to use'
        },
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

sub run_indexer {
    die("Must be defined in subclass");
}

sub prepare_gene_file {
    die("Must be defined in subclass");
}


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
        use_version => $self->picard_version,
        output_file => $sam_file,
        reference_fasta => $fasta_file,
        species => $species_name,
        genome_assembly => $assembly_name,
    )) {
        die('Failed to create sequence dictionary!');
    }

    # gmt picard sort-sam
    unless (Genome::Model::Tools::Picard::SortSam->execute(
        use_version => $self->picard_version,
        input_file => $sam_file,
        output_file => $seqdict_file,
    )) {
        die('Failed to sort sam file!');
    }
    # TODO:
    #die('Validate the order of the chromosome SQ lines with the reference build!');
    unlink($sam_file);
    return $seqdict_file;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return 'build_merged_alignments/chimerascan-index/' . $self->id;
}


1;
