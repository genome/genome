package Genome::InstrumentData::AlignmentResult::Command::CufflinksExpression;

use strict;
use warnings;

use Genome;
use Sys::Hostname;
use File::Path;

class Genome::InstrumentData::AlignmentResult::Command::CufflinksExpression {
    is => ['Command::V2'],
    has_input => [
        alignment_result_id => {
            is => 'Number',
            doc => 'ID of the result for the alignment data upon which to run coverage stats',
        },
        reference_build_id => {
            is => 'Number',
            doc => 'ID of the reference sequence used to generate metrics.',
        },
        annotation_build_id => {
            is => 'Number',
            doc => 'ID of the annotation used to generate metrics.',
            is_optional => 1,
        },
        expression_directory => {
            is => 'Text',
            doc => 'The directory to write cufflinks expression estimates.',
        },
    ],
    has_param => [
        cufflinks_params => {
            is => 'Text',
            doc => 'A string of extra params to pass to cufflinks.',
        },
        cufflinks_version => {
            is => 'Text',
            doc => 'The version of cufflinks to use.',
        },
        mask_reference_transcripts => {
            is => 'Text',
            doc => 'The name of an annotation file of known transcripts to mask during expression estimation.',
        },
        annotation_reference_transcripts_mode => {
            is => 'Text',
            doc => 'The mode for cufflinks to utilize the known annotation files.',
            valid_values => ['reference only', 'reference guided', 'de novo'],
        },
    ],
    has => [
        alignment_result => {
            is => 'Genome::SoftwareResult',
            id_by => 'alignment_result_id',
            doc => 'the alignment data upon which to run coverage stats',
        },
        reference_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            id_by => 'reference_build_id',
            doc => 'the reference sequence upon which to run',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            id_by => 'annotation_build_id',
            doc => 'the annotation upon which to run',
            is_optional => 1,
        },
    ],
};

sub execute {
    my $self = shift;

    # Reference inputs
    my $reference_build = $self->reference_build;
    my $reference_fasta_file = $reference_build->full_consensus_path('fa');
    unless (-s $reference_fasta_file) {
        $self->error_message("Reference FASTA File ($reference_fasta_file) is missing");
        return;
    }
    
    # Annotation inputs
    my $annotation_build = $self->annotation_build;

    # Alignment inputs
    my $alignment_result = $self->alignment_result;

    my $tophat_file;
    # Tophat v1.1.0 and later produces BAM output
    if (version->parse($alignment_result->aligner_version) >= version->parse('1.1.0')) {
        # Cufflinks v0.9.0 and later will work on sam or bam files
        if (version->parse($self->cufflinks_version) >= version->parse('0.9.0')) {
            $tophat_file = $alignment_result->bam_file;
        } else {
            $tophat_file = Genome::Sys->create_temp_file_path($self->build->id .'.sam');
            unless (Genome::Model::Tools::Sam::BamToSam->execute(
                bam_file => $alignment_result->bam_file,
                sam_file => $tophat_file,
            )) {
                $self->error_message('Failed to convert BAM '. $alignment_result->bam_file .' to tmp SAM file '. $tophat_file);
                return;
            }
        }
    } else {
        # Tophat versions before v1.1.0 and produce sam output and cufflinks __SHOULD__ run on sam or bam input
        $tophat_file = $alignment_result->sam_file;
    }
    
    my $expression_directory = $self->expression_directory;
    unless (-d $expression_directory) {
        Genome::Sys->create_directory($expression_directory);
    }

    my $params = $self->cufflinks_params || '';

    if (version->parse($self->cufflinks_version) >= version->parse('0.9.0')) {
        if (version->parse($self->cufflinks_version) < version->parse('1.0.0')) {
            $params .= ' -r '. $reference_fasta_file;
        } else {
            $params .= ' -b '. $reference_fasta_file;
        }
        if ($annotation_build) {
            # Determine the type of masking to use with Cufflinks
            my $mask_transcripts = $self->mask_reference_transcripts;
            if ($mask_transcripts) {
                my $mask_file_method = $mask_transcripts .'_file';
                my $mask_gtf_path = $annotation_build->$mask_file_method('gtf',$reference_build->id);
                unless(-s $mask_gtf_path) {
                    $mask_gtf_path = $annotation_build->$mask_file_method('gtf');
                }
                unless ($mask_gtf_path) {
                    $self->error_message('Failed to find GTF annotation used to mask transcripts with type: '. $mask_transcripts);
                    return;
                }
                $params .= ' -M '. $mask_gtf_path;
            }
            
            # Determine both the annotation file and mode to use it with Currlinks
            my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build->id);
            unless($gtf_path) {
                $gtf_path = $annotation_build->annotation_file('gtf');
            }
            my $mode = $self->annotation_reference_transcripts_mode;
            unless (defined($mode)) {
                $mode = 'de novo';
            }
            if ($mode eq 'reference only') {
                unless (defined($gtf_path)) {
                    $self->error_message('There is no annotation GTF file for annotation build: '. $annotation_build->id);
                    return;
                }
                $params .= ' -G '. $gtf_path;
            } elsif ($mode eq 'reference guided') {
                if (version->parse($self->cufflinks_version) < version->parse('1.0.0')) {
                    $self->error_message('This version of cufflinks '. $self->cufflinks_version .' is not compatible with the '. $mode .' annotation parameter -g.  Only v1.0.0 and greater!');
                    return;
                }
                unless (defined($gtf_path)) {
                    $self->error_message('There is no annotation GTF file found for annotation build: '. $annotation_build->id);
                    return;
                }
                $params .= ' -g '. $gtf_path;
            } elsif ($mode eq 'de novo') {
                # It doesn't matter if the files exist or not... we won't need it
            } else {
                $self->error_message('The processing_profile param annotation_reference_transcripts_mode \''. $mode .'\' is not supported!');
                return;
            }
        }
    }
    
    unless (Genome::Model::Tools::Cufflinks::Assemble->execute(
        input_file => $tophat_file,
        params => $params,
        output_directory => $expression_directory,
        use_version => $self->cufflinks_version,
    )) {
        $self->error_message('Failed to execute cufflinks!');
        return;
    }

    return 1;
}


1;
