package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels::Result;

use Genome;
use strict;
use warnings;
use Carp;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels::Result {
    is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase::Result',
    has_param => [
        varscan_version => { is => 'Text' },
    ],
};

sub _generate_result {
    my ($self, $staging_directory) = @_;
    my @builds = $self->builds;
    # FIXME pass command in, as non-input and non-param but required.
    my $cmd = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels->create(
            builds => \@builds,
            output_directory => $staging_directory,
            max_files_per_merge => $self->max_files_per_merge,
            roi_list => $self->roi_list,
            wingspan => $self->wingspan,
            allow_multiple_processing_profiles => $self->allow_multiple_processing_profiles,
            joinx_version => $self->joinx_version,
            varscan_version => $self->varscan_version,
    );
    my $return_value = $cmd->generate_result();
    Carp::croak($self->error_message('Command failed')) unless $return_value;
}
