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
        forced_variations_build_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'The ID of a Genome::Model::Build::ImportedVariationList containing a vcf of sites to force-genotype regardless of whether the model-group called it. ' .
                    'See - genome model imported-variation-list import-variants --help',
        },
    ],
    has => [
        variant_type => {
            is_constant => 1,
            value => 'snvs',
        }
    ],
    has_transient => [
        process => {
            is => 'Genome::Process',
            doc => 'The process wrapping the workflow that generates this result',
        },
    ],
};

sub execute {
    my ($self) = @_;

    Genome::Sys->create_directory($self->output_directory);
    $self->debug_message("Resolving Builds...");
    my $builds = $self->_resolve_builds();
    my %params = (
        forced_variations_build_id => $self->forced_variations_build_id,
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

    $params{command} = $self;
    $params{users} = {
        requestor => $self->_create_process(),
        sponsor   => Genome::Sys->current_user,
    };

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
        forced_variations_vcf => $self->_resolve_forced_variations_vcf,
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

# Make sure the file exists if provided. If not we need to fill in something so the workflow can proceed.
sub _resolve_forced_variations_vcf {
    my $self = shift;

    if ($self->forced_variations_build_id and $self->forced_variations_build_id ne $self->_undefined_forced_variations_build_id_value) {
        my $build = Genome::Model::Build::ImportedVariationList->get($self->forced_variations_build_id);
        unless ($build) {
            die $self->error_message("forced_variations_build_id provided (%s) fails to resolve a valid build", $self->forced_variations_build_id);
        }

        my $forced_variations_vcf = $build->snvs_vcf;
        Genome::Sys->validate_file_for_reading($forced_variations_vcf);
        return $forced_variations_vcf;
    } else {
        return $self->_undefined_forced_variations_build_id_value;
    }
}

sub _undefined_forced_variations_build_id_value {
        return 'NO_FORCED_VARIATIONS_BUILD';
}

1;
