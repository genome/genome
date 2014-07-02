package Genome::VariantReporting::Dbsnp::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Dbsnp::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        dbsnp_build_id => {
            is => 'Text',
        },
    ],
    has_param => [
        info_string => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
    ],
};


sub output_filename {
    return 'joinx_vcf_annotate.vcf.gz';
}

sub annotation_file {
    my $self = shift;

    my $dbsnp_build = Genome::Model::Build::ImportedVariationList->get($self->dbsnp_build_id);
    if ($dbsnp_build) {
        return $dbsnp_build->snvs_vcf;
    } else {
        die sprintf("No dbsnp build with id (%s)", $self->dbsnp_build_id);
    }
}

sub _run {
    my $self = shift;

    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);
    my $info_string = $self->info_string;
    my $info        = $info_string ? 1 : 0;

    my $vcf_annotator = Genome::Model::Tools::Joinx::VcfAnnotate->create(
        input_file      => $self->input_vcf,
        annotation_file => $self->annotation_file,
        output_file     => $output_file,
        use_bgzip       => 1,
        info_fields     => $info_string,
        info            => $info,
        use_version     => $self->joinx_version,
    );

    unless ($vcf_annotator->execute) {
        die $self->error_message("Failed to execute joinx vcf-annotate");
    }

    return 1;
}


1;
