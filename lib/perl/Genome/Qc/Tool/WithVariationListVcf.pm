package Genome::Qc::Tool::WithVariationListVcf;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::WithVariationListVcf {
    is => 'Genome::Qc::Tool',
    is_abstract => 1,
    has => [
        variation_list_build_id => {
            is => 'Text',
            is_optional => 1,
        },
        variation_list_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            id_by => 'variation_list_build_id',
        },
    ],
};

sub variation_list_vcf_file {
    my $self = shift;

    my $vcf_file = Genome::Sys->create_temp_file_path;
    my $cmd = Genome::Model::GenotypeMicroarray::Command::ExtractToVcf->create(
        variation_list_build => $self->variation_list_build,
        sample => $self->sample,
        output => $vcf_file,
    );
    unless ($cmd->execute) {
        die $self->error_message("ExtractToVcf command execution failed");
    }

    return $vcf_file;
}

1;
