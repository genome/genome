package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfSnvs;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfSnvs {
    is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase',
    has_optional_input => [
        samtools_version => {
            is => 'Text',
        },
        samtools_pileup_params => {
            is => 'Text',
            default => '-A -B',
        },
    ],
    has => [
        variant_type => {
            is_constant => 1,
            value => 'snvs',
        }
    ],
};

sub execute {
    my ($self) = @_;

    Genome::Sys->create_directory($self->output_directory);
    $self->debug_message("Resolving Builds...");
    my $builds = $self->_resolve_builds();
    my %params = (
        builds => $builds,
        max_files_per_merge => $self->max_files_per_merge,
        roi_list => $self->roi_list,
        wingspan => $self->wingspan,
        allow_multiple_processing_profiles => $self->allow_multiple_processing_profiles,
        joinx_version => $self->joinx_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
    unless($self->roi_list){
       delete $params{'roi_list'};
       delete $params{'wingspan'};
    }
    my $software_result = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfSnvs::Result->get_or_create(%params);
    $self->debug_message("Got or created software result with id "
        . $software_result->id . " (test_name='" . $software_result->test_name . "')");
    Genome::Sys->symlink_directory($software_result->output_dir,
        $self->output_directory);
    $self->software_result($software_result);
    $self->final_result(join("/", $software_result->output_dir,
                sprintf("%s.merged.vcf.gz", $self->variant_type)));
    return 1;
}

sub _get_workflow_xml {
    my $self = shift;
    my $xml_file;
    if($self->roi_list){
        $xml_file = "RegionLimitAndBackfillSnvVcf.xml";
    } else{
        $xml_file = "MergeAndBackfillSnvVcf.xml";
    }
    return $xml_file;
}

sub _get_variant_type_specific_inputs {
    my $self = shift;

    $self->_resolve_samtools_version();
    my %inputs = (
        samtools_version => $self->samtools_version,
        samtools_params => $self->samtools_pileup_params,
    );
    return \%inputs;
}

sub _resolve_samtools_version {
    my $self = shift;
    if (!defined($self->samtools_version)) {
        $self->samtools_version(Genome::Model::Tools::Sam->default_samtools_version);
    }
    return $self->samtools_version;
}

1;
