package Genome::Model::ReferenceSequence::Command::GenerateTranscriptome;

use warnings;
use strict;

use Genome;

class Genome::Model::ReferenceSequence::Command::GenerateTranscriptome{
    is => 'Genome::Command::Base',
    has_input => [
        reference_sequence_build_name => {
            is => 'Text',
            doc => 'the reference build, specified by name: ex. GRCh37-lite-build37',
        },
        annotation_build_name => {
            is => 'Text',
            doc => 'The annoation build, specified by name: ex. NCBI-human.ensembl/58_37c_v2',
        },
        merge_level => {
            doc => 'The transcriptome can consist of all transcript sequences or a squashed representation of each gene.',
            valid_values => ['gene','transcript'],
            default_value => 'gene',
        },
        output_fasta_file => {
            is => 'Text',
            doc => 'The output FASTA file that will contain all genes or transcript references in the forward strand orientation',
        },
        junctions_bed_file => {
            is => 'Text',
            doc => 'An output BED file that contains the splice junction positions(4bp) in the transcripts/genes.',
        },
    ],
    doc => 'Generate a transcriptome reference FASTA from a reference genome and annotation set(Ensembl-only for now).',
};

sub help_synopsis {
    my $class = shift;
    return <<EOS;
genome model reference-sequence generate-transcriptome --annotation-build-name='NCBI-human.ensembl/58_37c_v2' --reference-sequence-build-name=GRCh37-lite-build37 --merge-level=transcript --output-fasta-file=transcripts.fasta --junctions-bed-file=transcript_junctions.bed
EOS
}

sub help_detail {

    my $class = shift;
    return <<'EOS';
This tool is designed to take an annotation build and a reference sequence build and produce a transcriptome reference FASTA file.  Also a BED file of the location of splice junctions(2bp on either side) in context of the transcript/gene sequence coordinates.  All sequences are kept in a positive strand orientation relative to the genome reference.  The names of each sequence in the output FASTA include the gene or transcript id, the chromosome, and the strand in a pipe '|' separated string.
EOS
}

sub execute {
    my $self = shift;
    my $reference_sequence_build = Genome::Model::Build::ReferenceSequence->get(name => $self->reference_sequence_build_name);
    my $annotation_build = Genome::Model::Build::ImportedAnnotation->get(name => $self->annotation_build_name);
    my $bed_file;
    if ($self->merge_level eq 'gene') {
        # Get a squashed BED file with the strand
        $bed_file = $annotation_build->generate_annotation_file('bed',$reference_sequence_build->id,1,1);
    } elsif ($self->merge_level eq 'transcript') {
        # Get a complete BED file
        $bed_file = $annotation_build->generate_annotation_file('bed',$reference_sequence_build->id,0,0);
    } else {
        die('Something BAD!');
    }
    my $input_fasta_file = $reference_sequence_build->full_consensus_path('fa');
    unless (Genome::Model::Tools::RefCov::BuildReferenceFromRoi->execute(
        bed_file => $bed_file,
        stitch_level => $self->merge_level,
        junctions_bed_file => $self->junctions_bed_file,
        output_fasta_file => $self->output_fasta_file,
        input_fasta_file => $input_fasta_file,
        # The input BED file should serve as the coordinates of exons
        # coordinates => 1,
    )) {
        die('Failed to generate new FASTA!');
    }
    return 1;
}

1;
