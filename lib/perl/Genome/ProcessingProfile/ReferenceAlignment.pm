package Genome::ProcessingProfile::ReferenceAlignment;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::ReferenceAlignment {
    is => 'Genome::ProcessingProfile',
    has_param => [
        sequencing_platform => {
            doc => 'The sequencing platform from whence the model data was generated',
            valid_values => ['454', 'solexa', 'sanger'],
        },
        dna_type => {
            doc => 'the type of dna used in the reads for this model',
            valid_values => ['genomic dna', 'cdna']
        },
        transcript_variant_annotator_version => {
            doc => 'Version of the "annotate transcript-variants" tool to run during the annotation step',
            is_optional => 1,
            default_value => Genome::Model::Tools::Annotate::TranscriptVariants->default_annotator_version,
            valid_values => [ 0,1,2,3,4],#Genome::Model::Tools::Annotate::TranscriptVariants->available_versions ],
        },
        transcript_variant_annotator_filter => {
            doc => 'annotation-filter option to be used by the "annotate transcript-variants" tool run during the annotation step',
            is_optional => 1,
            default_value => 'top',
            valid_values => ['top', 'none', 'gene'],
        },
        transcript_variant_annotator_accept_reference_IUB_codes => {
            doc => 'annotation accept-reference-IUB-codes option to be used by the "annotate transcript-variants" to run during the annotation step',
            is_optional => 1,
            default_value => '0',
            valid_values => [0, 1],
        },
        alignment_strategy => {
            is => 'Text',
            is_many => 0,
            is_optional => 1,
            doc => 'Strategy to be used to align the instrument data (this will eventually replace many of the alignment parameters)',
        },
        snv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect snvs.",
        },
        indel_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect indels.",
        },
        sv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect svs.",
        },
        cnv_detection_strategy => {
            is => "Text",
            is_many => 0,
            is_optional =>1,
            doc => "Strategy to be used to detect cnvs.",
        },
        picard_version => {
            doc => 'picard version for MarkDuplicates, MergeSamfiles, CreateSequenceDictionary...',
            is_optional => 1,
        },
        samtools_version => {
            doc => 'samtools version for SamToBam, samtools merge, etc...',
            is_optional => 1,
        },
        bedtools_version => {
            doc => 'bedtools version for bedtools bamtofastq',
            is_optional => 1,
        },
        merger_name => {
            doc => 'name of bam merger, picard, samtools (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_version => {
            doc => 'version of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        merger_params => {
            doc => 'parameters of bam merger (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_name => {
            doc => 'name of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_version => {
            doc => 'version of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        duplication_handler_params => {
            doc => 'parameters of pcr duplication handler (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_name => {
            doc => 'alignment algorithm/software used for this model (this will be replaced by alignment_strategy)',
        },
        read_aligner_version => {
            doc => 'the aligner version used for this model (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        read_aligner_params => {
            doc => 'command line args for the aligner (this will be replaced by alignment_strategy)',
            is_optional => 1,
        },
        force_fragment => {
            is => 'Integer',
            #This doesn't seem to work yet because of the create code, can't the valid values logic be removed from create???
            default_value => '0',
            #valid_values => ['0', '1'],
            doc => 'force all alignments as fragment reads',
            is_optional => 1,
        },
        read_trimmer_name => {
            doc => 'trimmer algorithm/software used for this model',
            is_optional => 1,
        },
        read_trimmer_version => {
            doc => 'the trimmer version used for this model',
            is_optional => 1,
        },
        read_trimmer_params => {
            doc => 'command line args for the trimmer',
            is_optional => 1,
        },
        coverage_stats_params => {
            doc => 'parameters necessary for generating reference coverage in the form of two comma delimited lists split by a colon like 1,5,10,15,20:0,200,500',
            is_optional => 1,
        },
        append_event_steps => {
            doc => 'Event classes to append to event_stage_job_classes, e.g. "alignment => Genome::Model::Event::Build::ReferenceAlignment::QC::CopyNumber".',
            is_optional => 1,
            is_deprecated => 1,
        },

    ],
};

sub _resolve_type_name_for_class {
    return 'reference alignment';
}


1;
