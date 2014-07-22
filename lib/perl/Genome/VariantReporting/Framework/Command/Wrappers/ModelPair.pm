package Genome::VariantReporting::Framework::Command::Wrappers::ModelPair;

use strict;
use warnings;

use Genome;

use File::Spec;

class Genome::VariantReporting::Framework::Command::Wrappers::ModelPair {
    has => {
        discovery => { is => 'Genome::Model::Build', },
        validation => { is => 'Genome::Model::Build', },
        base_output_dir => { is => 'Text', },
    },
    has_calculated => {
        output_dir => {
            calculate_from => [qw/ base_output_dir discovery /],
            calculate => q| my $roi_nospace  = $discovery->region_of_interest_set_name;
                $roi_nospace =~ s/ /_/g;
                return File::Spec->join($base_output_dir, $roi_nospace); |,
        },
        resource_file => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "resource.yaml") ),
        },
    },
};

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
    $self->create_input_vcf("snvs");
    $self->create_input_vcf("indels");
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

sub generate_resource_file {
    my $self = shift;

    return if not $self->is_valid;
    my $resource = {};

    my @aligned_bams;
    push @aligned_bams, $self->discovery->merged_alignment_result->id;
    push @aligned_bams, $self->discovery->control_merged_alignment_result->id;
    push @aligned_bams, $self->validation->control_merged_alignment_result->id;
    $resource->{aligned_bam_result_id} = \@aligned_bams;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $self->discovery->region_of_interest_set_name);
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

sub input_vcf {
    my ($self, $variant_type) = @_;
    return File::Spec->join($self->output_dir, "$variant_type.vcf.gz");
}

sub create_input_vcf {
    my ($self, $variant_type) = @_;
    my $vcf = $self->discovery->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
    my $sed_cmd = sprintf("zcat %s | %s | gzip > %s", $vcf, 'sed s/^#CHROM.*/\&\	'.$self->validation->tumor_sample->name."/", $self->input_vcf($variant_type));
    my $rv = Genome::Sys->shellcmd(
        cmd => $sed_cmd,
        output_files => [$self->input_vcf($variant_type)],
        input_files => [$vcf],
    );
}
1;

