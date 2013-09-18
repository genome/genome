package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::RegionLimitVcf;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::CreateCrossSampleVcf::RegionLimitVcf {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build',
        },
        variant_type => {
            is => 'Text',
            default => 'snvs',
            valid_values => ['snvs','indels'],
            doc => 'The type of variations present in the vcf files.',
        },
        output_directory => {
            is => 'Path',
        },
        region_bed_file => {
            is => 'Text',
            doc => 'The bed coordinates of regions to limit the vcf to',
        },
        roi_name => {
            is => 'Text',
            doc => 'Region of Interest Set Name, this will be added as a tag in the vcf header is specified.',
            is_optional => 1,
        },
        wingspan => {
            is => 'Text',
            doc => 'This is the amount of bases to include on either side of each region',
            default => 0,
        },
    ],
    has_optional_output => [
        output_vcf => {
            is => 'Path',
        },
    ],
};

sub execute {
    my $self = shift;

    $self->output_vcf($self->output_file);
    my $cmd = Genome::Model::Tools::Vcf::RegionLimit->create(
        output_file => $self->output_vcf,
        vcf_file => $self->vcf_file,
        region_bed_file => $self->region_bed_file,
        roi_name => $self->roi_name,
        wingspan => $self->wingspan,
    );
    return $cmd->execute();
}

sub output_file {
    my $self = shift;

    my $filename = sprintf("%s.%s.region_limited.vcf.gz",
        $self->variant_type, $self->build->model->subject->id);
    return File::Spec->join($self->output_directory, $filename);
}

sub vcf_file {
    my $self = shift;

    my $accessor = sprintf("get_%s_vcf", $self->variant_type);
    return $self->build->$accessor;
}

1;
