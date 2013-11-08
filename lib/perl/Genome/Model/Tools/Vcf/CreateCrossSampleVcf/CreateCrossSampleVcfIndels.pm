package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels {
    is => 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfBase',
    has => [
        variant_type => {
            is_constant => 1,
            value => 'indels',
        },
        varscan_version => {
            is => 'Text',
            doc => 'Varscan version to use in all varscan operations',
            default => '2.3.6', # TODO lean on gmt varscan for its default (same with joinx version)
        },
    ],
};

sub execute {
    my ($self) = @_;

    $self->status_message("Resolving Builds...");
    my $builds = $self->_resolve_builds();
    my %params = (
        builds => $builds,
        max_files_per_merge => $self->max_files_per_merge,
        roi_list => $self->roi_list,
        wingspan => $self->wingspan,
        allow_multiple_processing_profiles => $self->allow_multiple_processing_profiles,
        joinx_version => $self->joinx_version,
        varscan_version => $self->varscan_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
    );
    unless($self->roi_list){
       delete $params{'roi_list'};
       delete $params{'wingspan'};
    }
    my $software_result = Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcfIndels::Result->get_or_create(%params);
    $self->status_message("Got or created software result with id "
        . $software_result->id . " (test_name='" . $software_result->test_name . "')");
    Genome::Sys->symlink_directory($software_result->output_dir,
        $self->output_directory);
    $self->software_result($software_result);
    $self->final_result(join("/", $software_result->output_dir,
                sprintf("%s.merged.vcf.gz", $self->variant_type)));
    return 1;
}


sub _get_workflow_inputs {
    my $self = shift;
    my ($builds, $variant_type_specific_inputs, $region_limiting_specific_inputs) = @_;
    my $inputs = $self->SUPER::_get_workflow_inputs(@_);
    $inputs->{varscan_version} = $self->varscan_version;
    return $inputs;
}

sub _get_workflow_xml {
    my $self = shift;
    my $xml_file;
    if ($self->roi_list){
        $xml_file = "RegionLimitAndBackfillIndelVcf.xml";
    } else{
        $xml_file = "MergeAndBackfillIndelVcf.xml";
    }
    return $xml_file;
}

sub _get_variant_type_specific_inputs {
    my $self = shift;

    my @input_bams = map {$_->whole_rmdup_bam_file} $self->builds;
    my %inputs = (
        input_bams => \@input_bams,
        output_directory => File::Spec->join($self->output_directory, 'indel_backfilling'),
        exact_pos => 1,
    );
    return \%inputs;
}

1;
