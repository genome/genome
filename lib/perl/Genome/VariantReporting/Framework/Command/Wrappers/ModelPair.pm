package Genome::VariantReporting::Framework::Command::Wrappers::ModelPair;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::VariantReporting::Framework::Command::Wrappers::ModelPair {
    has => {
        discovery => { is => 'Genome::Model::Build', },
        validation => { is => 'Genome::Model::Build', },
        roi => { is => 'Text', },
        base_output_dir => { is => 'Text', },
    },
    has_calculated => {
        output_dir => {
            calculate_from => [qw/ base_output_dir roi /],
            calculate => q| return File::Spec->join($base_output_dir, $roi); |,
        },
        resource_file => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "resource.yaml") ),
        },
        reports_directory => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "reports"); ),
        },
        logs_directory => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "logs"); ),
        },
    },
};

sub is_valid {
    my $self = shift;

    if (my @problems = $self->__errors__) {
        $self->error_message('Model pair is invalid!');
        for my $problem (@problems) {
            my @properties = $problem->properties;
            $self->error_message("Property " .
                join(',', map { "'$_'" } @properties) .
                ': ' . $problem->desc);
        }
        return;
    }

    return 1;
}

sub generate_resource_file {
    my $self = shift;

    return if not $self->is_valid;
    my $resource = {};

    my @aligned_bams;
    push @aligned_bams, $self->discovery->merged_alignment_result->id;
    push @aligned_bams, $self->discovery->control_merged_alignment_result->id;
    push @aligned_bams, $self->validation->control_merged_alignment_result->id;
    $resource->{aligned_bams} = \@aligned_bams;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $self->roi);
    $feature_list_ids{ON_TARGET} = $on_target_feature_list->id;
    my $segdup_feature_list = $self->discovery->get_feature_list("segmental_duplications");
    $feature_list_ids{SEG_DUP} = $segdup_feature_list->id;
    $resource->{feature_list_ids} = \%feature_list_ids;

    $resource->{reference_fasta} = $self->discovery->reference_sequence_build->full_consensus_path("fa");

    my %translations;
    $translations{d0_tumor} = $self->discovery->tumor_sample->name;
    $translations{d30_normal} = $self->discovery->normal_sample->name;
    $translations{d30_tumor} = $self->validation->tumor_sample->name;
    $resource->{translations} = \%translations;

    $resource->{dbsnp_vcf} = $self->discovery->previously_discovered_variations_build->snvs_vcf;

    YAML::DumpFile(File::Spec->join($self->resource_file), $resource);

    return 1;
}

1;

