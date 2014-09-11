package Genome::VariantReporting::Framework::Command::Wrappers::ModelPair;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;

class Genome::VariantReporting::Framework::Command::Wrappers::ModelPair {
    has => {
        discovery => { is => 'Genome::Model::Build', },
        validation => {
            is => 'Genome::Model::Build',
            is_optional => 1,
        },
        base_output_dir => { is => 'Text', },
        plan_file_basename => {
            is => 'Text',
            default_value => "cle_TYPE_report.yaml",
            doc => "plan file name where 'snvs' or 'indels' is substituted by placeholder TYPE",
        },
    },
    has_calculated => {
        output_dir => {
            calculate_from => [qw/ base_output_dir discovery /],
            calculate => q| my $roi_nospace  = $discovery->region_of_interest_set->name;
                $roi_nospace =~ s/ /_/g;
                return File::Spec->join($base_output_dir, $roi_nospace); |,
        },
        resource_file => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "resource.yaml") ),
        },
    },
};

sub plan_file {
    my ($self, $type) = @_;
    my $base_name = $self->plan_file_basename;
    $base_name =~ s/TYPE/$type/;
    return File::Spec->join($self->_plan_search_dir, $base_name);
}

sub _plan_search_dir {
    my $variant_reporting_base_dir = dirname(dirname(dirname(dirname(__FILE__))));
    return File::Spec->join($variant_reporting_base_dir, 'plan_files');
}

sub report_names {
    return qw(cle_full_report cle_simple_report);
}

sub reports_directory {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_dir, "reports_$variant_type");
};

sub logs_directory {
    my ($self, $variant_type) = @_;
    return  File::Spec->join($self->output_dir, "logs_$variant_type");
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    Genome::Sys->create_directory($self->output_dir);
    Genome::Sys->create_directory($self->reports_directory("snvs"));
    Genome::Sys->create_directory($self->reports_directory("indels"));
    Genome::Sys->create_directory($self->logs_directory("snvs"));
    Genome::Sys->create_directory($self->logs_directory("indels"));
    $self->generate_resource_file;
    return $self;
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

sub get_aligned_bams {
    my $self = shift;
    my @aligned_bams;
    push @aligned_bams, $self->discovery->merged_alignment_result->id;
    push @aligned_bams, $self->discovery->control_merged_alignment_result->id;
    if ($self->validation) {
        push @aligned_bams, $self->validation->merged_alignment_result->id;
    }
    return \@aligned_bams;
}

sub get_translations {
    my $self = shift;
    my %translations;
    $translations{d0_tumor} = $self->discovery->tumor_sample->name;
    $translations{d30_normal} = $self->discovery->normal_sample->name;
    if ($self->validation) {
        $translations{d30_tumor} = $self->validation->tumor_sample->name;
    }
    return \%translations;
}

sub generate_resource_file {
    my $self = shift;

    return if not $self->is_valid;
    my $resource = {};

    $resource->{aligned_bam_result_id} = $self->get_aligned_bams;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $self->discovery->region_of_interest_set->name);
    $feature_list_ids{ON_TARGET} = $on_target_feature_list->id;
    my $segdup_feature_list = $self->discovery->get_feature_list("segmental_duplications");
    $feature_list_ids{SEG_DUP} = $segdup_feature_list->id;
    # TODO: There has to be a better way...
    $feature_list_ids{AML_RMG} = '0e4973c600244c3f804d54bee6f81145';
    $resource->{feature_list_ids} = \%feature_list_ids;

    $resource->{reference_fasta} = $self->discovery->reference_sequence_build->full_consensus_path("fa");

    $resource->{translations} = $self->get_translations;

    $resource->{dbsnp_vcf} = $self->discovery->previously_discovered_variations_build->snvs_vcf;

    YAML::DumpFile(File::Spec->join($self->resource_file), $resource);

    return 1;
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    return $self->discovery->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
}
1;

