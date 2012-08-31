package Genome::Model::Tools::Vcf::CreateCrossSampleVcfResult;

use Genome;
use strict;
use warnings;
use Carp;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcfResult {
    is => 'Genome::SoftwareResult::DiskAllocationStaged',

    has_input => [
        builds => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            is_many => 1,
        },
    ],
    has_param => [
        max_files_per_merge => { is => 'Text' },
        variant_type => { is => 'Text' },
        roi_file => { is => 'Text' },
        roi_name => { is => 'Text' },
        wingspan => { is => 'Text' },
        allow_multiple_processing_profiles => { is => 'Boolean' },
        joinx_version => { is => 'Text' },
    ],
};

sub _generate_result {
    my ($self, $staging_directory) = @_;
    my @builds = $self->builds;
    my $cmd = Genome::Model::Tools::Vcf::CreateCrossSampleVcf->create(
            builds => \@builds,
            output_directory => $staging_directory,
            max_files_per_merge => $self->max_files_per_merge,
            variant_type => $self->variant_type,
            roi_file => $self->roi_file,
            roi_name => $self->roi_name,
            wingspan => $self->wingspan,
            allow_multiple_processing_profiles => $self->allow_multiple_processing_profiles,
            joinx_version => $self->joinx_version,
    );
    Carp::croak($self->error_message('Command failed')) unless $cmd->generate_result();
}

1;
