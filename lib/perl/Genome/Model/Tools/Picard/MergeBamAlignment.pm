package Genome::Model::Tools::Picard::MergeBamAlignment;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Picard::MergeBamAlignment {
    is  => 'Genome::Model::Tools::Picard',
    has_input => [
        unmapped_bam => {
            is => 'Text',
            doc => 'Original SAM or BAM file of unmapped reads, which must be in queryname order.',
        },
        aligned_bam => {
            is => 'Text',
            doc => 'SAM or BAM file with alignment data',
        },
        output_file => {
            is => 'Text',
            doc => 'Merged SAM or BAM file to write to.',
        },
        reference_sequence => {
            is => 'Text',
            doc => 'Path to the fasta file for the reference sequence.',
            example_values => ['/gscmnt/gc4096/info/model_data/2741951221/build101947881/all_sequences.fa'],
        },
        paired_run => {
            is => 'Boolean',
            doc => 'Whether this is a paired-end run.',
            default_value => 1,
            is_optional => 1,
        },
        max_insertions_or_deletions => {
            is => 'Integer',
            doc => 'The maximum number of insertions or deletions permitted for an alignment to be included. Alignments with more than this many insertions or deletions will be ignored. Default value: 1. This option can be set to \'null\' to clear the default value.',
            default_value => '0',
            is_optional => 1,
        },
        primary_alignment_strategy =>{
            is => 'Text',
            doc => 'Strategy for selecting primary alignment when the aligner has provided more than one alignment for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary alignment is filtered out for some reason. BestMapq expects that multiple alignments will be correlated with HI tag, and prefers the pair of alignments with the largest MAPQ, in the absence of a primary selected by the aligner. EarliestFragment prefers the alignment which maps the earliest base in the read. Note that EarliestFragment may not be used for paired reads. BestEndMapq is appropriate for cases in which the aligner is not pair-aware, and does not output the HI tag. It simply picks the alignment for each end with the highest MAPQ, and makes those alignments primary, regardless of whether the two alignments make sense together.MostDistant is also for a non-pair-aware aligner, and picks the alignment pair with the largest insert size. If all alignments would be chimeric, it picks the alignments for each end with the best MAPQ. For all algorithms, ties are resolved arbitrarily. Default value: BestMapq.',
            valid_values => ['BestMapq', 'EarliestFragment', 'BestEndMapq', 'MostDistant'],
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'Tool to merge BAM or SAM alignment files using Picard';
}

sub help_detail {
    return <<EOS
     Merges alignment data from a SAM or BAM file with additional data stored in an unmapped BAM file and produces a third SAM or BAM file of aligned and unaligned reads. NOTE that this program expects to find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it to have the same base name as the reference fasta except with the extension '.dict'

     For Picard documentation of this command see:
     http://picard.sourceforge.net/command-line-overview.shtml#MergeBamAlignment
EOS
}

sub execute {
    my $self = shift;
    if ($self->use_version < 1.23) {
        die('Please use Picard version 1.23 or greater to run '. __PACKAGE__);
    }
    my $merge_cmd = $self->picard_path .'/MergeBamAlignment.jar net.sf.picard.sam.MergeBamAlignment  UNMAPPED_BAM='.
        $self->unmapped_bam .' ALIGNED_BAM='. $self->aligned_bam .' OUTPUT='. $self->output_file
        .' REFERENCE_SEQUENCE='. $self->reference_sequence .' MAX_INSERTIONS_OR_DELETIONS='. $self->max_insertions_or_deletions;
    if ($self->paired_run) {
        $merge_cmd .= ' PAIRED_RUN=true';
    } else {
        $merge_cmd .= ' PAIRED_RUN=false';
    }
    if (defined($self->primary_alignment_strategy)) {
        $merge_cmd .= ' PRIMARY_ALIGNMENT_STRATEGY='. $self->primary_alignment_strategy;
    }
    $self->run_java_vm(
        cmd => $merge_cmd,
        input_files => [$self->unmapped_bam, $self->aligned_bam, $self->reference_sequence],
        output_files => [$self->output_file],
        skip_if_output_is_present => 0,
    );
    return 1;
}


1;
