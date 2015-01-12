package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase::Result;

use Genome;
use strict;
use warnings;
use Carp;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase::Result {
    is => 'Genome::SoftwareResult::DiskAllocationStaged',
    is_abstract => 1,

    has_input => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
        },
    ],
    has_param => [
        max_files_per_merge => { is => 'Text' },
        roi_list => { is => 'Genome::FeatureList', is_optional => 1 },
        forced_variations_build_id => { is => 'Text', is_optional => 1 },
        wingspan => { is => 'Text', is_optional => 1 },
        allow_multiple_processing_profiles => { is => 'Boolean' },
        joinx_version => { is => 'Text' },
    ],
    has_transient => [
        command => {
            is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase',
            doc => 'The command for constructing this result',
        },
    ],
};

sub _generate_result {
    my ($self, $staging_directory) = @_;
    my @builds = $self->builds;

    my $cmd = $self->command or die ('No command');

    my $return_value = $cmd->generate_result($staging_directory);
    Carp::croak($self->error_message('Command failed')) unless $return_value;
}

1;
