package Genome::Model::ClinSeq::Command::AnnotateVcf::Result;

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);

class Genome::Model::ClinSeq::Command::AnnotateVcf::Result {
    is => 'Genome::SoftwareResult::StageableSimple::SingleFile',
    has_input => [
        input_file => {
            is => 'Text',
            doc => 'Vcf File to filter',
        },
        annotation_file => {
            is => 'Text',
            doc => 'Vcf File containing annotation',
        },
    ],
    has_optional_input => [
        info_fields => {
            is => 'Text',
            doc => 'Field ids to embed from the annotation VCF. Use colons to separate multiple field descriptors.',
            #doing the above because UR autosplits on commas with is_many, but joinx uses commas in its field descriptors
        },
        identifiers => {
            is => 'Boolean',
            doc => 'copy identifiers from the annotation file',
        },
        info => {
            is => 'Boolean',
            doc => 'copy information from info fields from the annotation file',
        },
    ],
};

sub _run {
    my $self = shift;

    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file => $self->input_file,
        annotation_file => $self->annotation_file,
        output_file => $self->_temp_staging_file_path,
        use_bgzip => $self->use_bgzip,
        info_fields => $self->info_fields,
        info => $self->info,
        use_version => Genome::Model::Tools::Joinx->get_default_version,
    );
    unless ($vcf_annotator->execute) {
        $self->fatal_message("Failed to execute Joinx Vcf annotation for vcf file (%s) using db (%s)", $self->input_file, $self->annotation_file);
    }

    return 1;
}

sub use_bgzip {
    my $self = shift;
    if ($self->input_file =~ /\.gz$/) {
        return 1;
    }
    else {
        return 0;
    }
}

sub _file_name {
    my $self = shift;
    return basename($self->input_file);
}

1;
