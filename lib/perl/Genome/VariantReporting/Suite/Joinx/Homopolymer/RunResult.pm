package Genome::VariantReporting::Suite::Joinx::Homopolymer::RunResult;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;

class Genome::VariantReporting::Suite::Joinx::Homopolymer::RunResult {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Result',
    has_input => [
        homopolymer_list_id => {
            is => 'String',
        },
    ],
    has_param => [
        info_string => {
            is => 'Text',
        },
        joinx_version => {
            is => 'Text',
        },
        max_length => {
            is => 'Integer',
        },
    ],
};


sub output_filename {
    return 'joinx_vcf_annotate_homopolymer.vcf.gz';
}

sub _run {
    my $self = shift;
    
    my $list_id          = $self->homopolymer_list_id;
    my $feature_list     = Genome::FeatureList->get($list_id);
    my $homopolymer_file = $feature_list->get_tabix_and_gzipped_bed_file;

    my $output_file = File::Spec->join($self->temp_staging_directory, $self->output_filename);

    my %params = (
        input_file       => $self->input_vcf,
        info_field       => $self->info_string,
        max_length       => $self->max_length,  
        use_version      => $self->joinx_version,
        output_file      => $output_file,
        use_bgzip        => 1,
        homopolymer_file => $homopolymer_file,
    );

    my $homopolymer_annotator = Genome::Model::Tools::Joinx::VcfAnnotateHomopolymer->create(%params);
    unless ($homopolymer_annotator->execute) {
        die $self->error_message("Failed to execute joinx vcf-annotate-homopolymer");
    }

    return 1;
}


1;
