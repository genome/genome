package Genome::VariantReporting::Command::Wrappers::ModelPair;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;
use Genome::File::Tsv;

class Genome::VariantReporting::Command::Wrappers::ModelPair {
    has => {
        discovery => { is => 'Genome::Model::Build', },
        followup => {
            is => 'Genome::Model::Build',
            is_optional => 1,
        },
        gold_sample_name => {
            is => 'Text',
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
        translations_file => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "resource.yaml") ),
        },
        sample_legend => {
            calculate_from => [qw/ output_dir /],
            calculate => q( File::Spec->join($output_dir, "sample_legend.tsv") ),
        }
    },
};

sub plan_file {
    my ($self, $type) = @_;
    my $base_name = $self->plan_file_basename;
    $base_name =~ s/TYPE/$type/;
    return File::Spec->join($self->_plan_search_dir, $base_name);
}

sub _plan_search_dir {
    my $variant_reporting_base_dir = dirname(dirname(dirname(__FILE__)));
    return File::Spec->join($variant_reporting_base_dir, 'plan_files');
}

sub report_names {
    my ($self, $variant_type) = @_;
    my $plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_file($self->plan_file($variant_type));
    my @file_names;
    for my $reporter_plan ($plan->reporter_plans) {
        push @file_names, $reporter_plan->object->file_name;
    }
    return @file_names;
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
    $self->generate_sample_legend_file;
    $self->generate_translations_file;
    my $provider = Genome::VariantReporting::Framework::Component::RuntimeTranslations->create_from_file($self->translations_file);
    for my $variant_type (qw(snvs indels)) {
        Genome::Sys->create_directory($self->reports_directory($variant_type));
        my $plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_file($self->plan_file($variant_type));
        $plan->write_to_file(File::Spec->join($self->reports_directory($variant_type), "plan.yaml"));
        Genome::Sys->create_directory($self->logs_directory($variant_type));
        $provider->write_to_file(File::Spec->join($self->reports_directory($variant_type), "resources.yaml"));
    }
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
    if ($self->followup) {
        push @aligned_bams, $self->followup->merged_alignment_result->id;
    }
    return \@aligned_bams;
}

sub get_sample_and_bam_map {
    my $self = shift;

    my %bams = (
        $self->discovery->tumor_sample->name  => $self->discovery->tumor_bam,
        $self->discovery->normal_sample->name => $self->discovery->normal_bam,
    );
    if ($self->followup) {
        $bams{$self->followup->tumor_sample->name} = $self->followup->tumor_bam,
    }
    return %bams;
}

sub get_translations {
    my $self = shift;
    my %translations;
    $translations{discovery_tumor} = $self->discovery->tumor_sample->name;
    $translations{normal} = $self->discovery->normal_sample->name;
    if ($self->followup) {
        $translations{followup_tumor} = $self->followup->tumor_sample->name;
    }
    if ($self->gold_sample_name) {
        $translations{gold} = $self->gold_sample_name;
    }
    return \%translations;
}

sub generate_sample_legend_file {
    my $self = shift;

    my $translations = $self->get_translations;
    my $legend_file = Genome::File::Tsv->create($self->sample_legend);
    my @headers = ('sample label', 'sample name');
    my $writer = $legend_file->create_writer(headers => \@headers);
    while ( my ($sample_label, $sample_name) = each %$translations ) {
        $writer->write_one({
            'sample label' => $sample_label,
            'sample name' => $sample_name,
         });
    }
}

sub generate_translations_file {
    my $self = shift;

    return if not $self->is_valid;

    my $translations = $self->get_translations;

    $translations->{aligned_bam_result_id} = $self->get_aligned_bams;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $self->discovery->region_of_interest_set->name);
    $feature_list_ids{ON_TARGET} = $on_target_feature_list->id;
    my $segdup_feature_list = $self->discovery->get_feature_list("segmental_duplications");
    $feature_list_ids{SEG_DUP} = $segdup_feature_list->id;
    # TODO: There has to be a better way...
    $feature_list_ids{AML_RMG} = '0e4973c600244c3f804d54bee6f81145';
    $translations->{feature_list_ids} = \%feature_list_ids;
    $translations->{homopolymer_list_id} = '696318bab30d47d49fab9afa845691b7';

    $translations->{reference_fasta} = $self->reference_sequence_build->full_consensus_path("fa");

    $translations->{dbsnp_vcf} = $self->discovery->previously_discovered_variations_build->snvs_vcf;
    $translations->{nhlbi_vcf} = _get_nhlbi_vcf(); 

    YAML::DumpFile(File::Spec->join($self->translations_file), $translations);

    return 1;
}

sub reference_sequence_build {
    my $self = shift;
    return $self->discovery->reference_sequence_build;
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    return $self->discovery->get_detailed_vcf_result($variant_type)->get_vcf($variant_type);
}

sub _get_nhlbi_vcf {
    return Genome::Model::Build::ImportedVariationList->get(
        version    => '2012.07.23', 
        model_name => 'nhlbi-esp-GRCh37-lite-build37',
    )->snvs_vcf;
}

1;

