package Genome::Model::Event::Build::RnaSeq::Expression::Cufflinks;

use strict;
use warnings;

use version;

use Genome;

class Genome::Model::Event::Build::RnaSeq::Expression::Cufflinks {
    is => ['Genome::Model::Event::Build::RnaSeq::Expression'],
    has => [
    ],
};

sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64 && mem>=32000] rusage[mem=32000] span[hosts=1]' -M 32000000 -n 4";
}

sub execute {
    my $self = shift;
    my $expression_directory = $self->build->accumulated_expression_directory;
    unless (-d $expression_directory) {
        Genome::Sys->create_directory($expression_directory);
    }
    my $alignment_result = $self->build->alignment_result;
    my $tophat_file;
    # Tophat v1.1.0 and later produces BAM output
    if (version->parse($alignment_result->aligner_version) >= version->parse('1.1.0')) {
        # Cufflinks v0.9.0 and later will work on sam or bam files
        if (version->parse($self->model->expression_version) >= version->parse('0.9.0')) {
            $tophat_file = $alignment_result->bam_file;
        } else {
            $tophat_file = Genome::Sys->create_temp_file_path($self->build->id .'.sam');
            unless (Genome::Model::Tools::Sam::BamToSam->execute(
                bam_file => $alignment_result->bam_file,
                sam_file => $tophat_file,
            )) {
                $self->error_message('Failed to convert BAM '. $alignment_result->bam_file .' to tmp SAM file '. $tophat_file);
                die($self->error_message);
            }
        }
    } else {
        # Tophat versions before v1.1.0 and produce sam output and cufflinks __SHOULD__ run on sam or bam input
        $tophat_file = $alignment_result->sam_file;
    }
    my $params = $self->model->expression_params || '';
    if (version->parse($self->model->expression_version) >= version->parse('0.9.0')) {
        my $reference_build = $self->model->reference_sequence_build;
        my $reference_path = $reference_build->full_consensus_path('fa');
        if (version->parse($self->model->expression_version) < version->parse('1.0.0')) {
            $params .= ' -r '. $reference_path;
        } else {
            $params .= ' -b '. $reference_path;
        }
        my $annotation_reference_transcripts = $self->model->annotation_reference_transcripts;
        if ($annotation_reference_transcripts) {
            my ($annotation_name,$annotation_version) = split(/\//, $annotation_reference_transcripts);
            my $annotation_model = Genome::Model->get(name => $annotation_name);
            unless ($annotation_model){
                $self->error_message('Failed to get annotation model for annotation_reference_transcripts: ' . $annotation_reference_transcripts);
                return;
            }

            unless (defined $annotation_version) {
                $self->error_message('Failed to get annotation version from annotation_reference_transcripts: '. $annotation_reference_transcripts);
                return;
            }

            my $annotation_build = $annotation_model->build_by_version($annotation_version);
            unless ($annotation_build){
                $self->error_message('Failed to get annotation build from annotation_reference_transcripts: '. $annotation_reference_transcripts);
                return;
            }

            # Determine the type of masking to use with Cufflinks
            my $mask_transcripts = $self->model->mask_reference_transcripts;
            if ($mask_transcripts) {
                my $mask_file_method = $mask_transcripts .'_file';
                my $mask_gtf_path = $annotation_build->$mask_file_method('gtf',$reference_build->id);
                unless ($mask_gtf_path) {
                    die('Failed to find GTF annotation used to mask transcripts with type: '. $mask_transcripts);
                }
                $params .= ' -M '. $mask_gtf_path;
            }

            # Determine both the annotation file and mode to use it with Currlinks
            my $gtf_path = $annotation_build->annotation_file('gtf',$reference_build->id);
            my $mode = $self->model->annotation_reference_transcripts_mode;
            unless (defined($mode)) {
                $mode = 'de novo';
            }
            if ($mode eq 'reference only') {
                unless (defined($gtf_path)) {
                    die('There is no annotation GTF file defined for annotation_reference_transcripts build: '. $annotation_reference_transcripts);
                }
                $params .= ' -G '. $gtf_path;
            } elsif ($mode eq 'reference guided') {
                if (version->parse($self->model->expression_version) < version->parse('1.0.0')) {
                    die('This version of cufflinks '. $self->model->expression_version .' is not compatible with the '. $mode .' annotation parameter -g.  Only v1.0.0 and greater!');
                }
                unless (defined($gtf_path)) {
                    die('There is no annotation GTF file defined for annotation_reference_transcripts build: '. $annotation_reference_transcripts);
                }
                $params .= ' -g '. $gtf_path;
            } elsif ($mode eq 'de novo') {
                # It doesn't matter if the files exist or not... we won't need it
            } else {
                die('The processing_profile param annotation_reference_transcripts_mode \''. $mode .'\' is not supported!');
            }
        }
        # Cufflinks should probably run once with the annotation gtf and -G to identify known transcripts.
        # Cufflinks should also run a second time to identify novel transcripts
        # This could be a param in the processing profile; however, resolving the annotation set(hence the right build/version) is performed here...
        # jwalker 06/24/2011 - You can easily define any set of params in the processing profile and not define annotation_reference_transcripts for the model
    }
    unless (Genome::Model::Tools::Cufflinks::Assemble->execute(
        input_file => $tophat_file,
        params => $params,
        output_directory => $expression_directory,
        use_version => $self->model->expression_version,
    )) {
        $self->error_message('Failed to execute cufflinks!');
        die($self->error_message);
    }
    return 1;
}

sub verify_successful_completion {
    my $self = shift;
    warn ('Please implement vsc for class '. __PACKAGE__);
    return 1;
}

1;
