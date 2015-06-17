package Genome::Qc::Tool::WithVariationListVcf;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::WithVariationListVcf {
    is_abstract => 1,
    has => [
        variation_list_build_id => {
            is => 'Text',
            is_optional => 1,
        }
    ],
};

sub variation_list_vcf_file {
    my $self = shift;
    my $variation_list_build = Genome::Model::Build::ImportedVariationList->get($self->variation_list_build_list);
    my $vcf_file = Genome::Sys->create_temp_file_path;
    my $cmd = Genome::Model::GenotypeMicroarray::Command::ExtractToVcf->create(
        variation_list_build => $variation_list_build,
        sample => $self->sample,
        output => $vcf_file,
    );
    unless ($cmd->execute) {
        die $self->error_message("ExtractToVcf command execution failed");
    }
    return $vcf_file;
}

1;
