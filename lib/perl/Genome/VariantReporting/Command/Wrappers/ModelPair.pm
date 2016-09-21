package Genome::VariantReporting::Command::Wrappers::ModelPair;

use strict;
use warnings;

use Genome;

use File::Basename qw(dirname);
use File::Spec;
use Genome::File::Tsv;
use YAML qw();
use List::MoreUtils qw(uniq);

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
        label => { is => 'Text', },
        plan_file_basename => {
            is => 'Text',
            default_value => "cle_TYPE_report.yaml",
            doc => "plan file name where 'snvs' or 'indels' is substituted by placeholder TYPE",
        },
    },
    has_transient => [
        translations_file => {
            is => 'Path',
        }
    ],
    has_transient_optional => {
        dag => {
            is => 'Genome::WorkflowBuilder::DAG',
        },
        common_translations => {
            is => 'HASH',
            default => {},
        },
    }
};

sub dag {
    my $self = shift;

    unless (defined($self->__dag)) {
        my $cmd = Genome::VariantReporting::Command::CreateMergedReports->create(
            %{$self->params_for_command},
        );
        my $dag = $cmd->dag;
        $dag->name(sprintf('%s (%s-%s)',
                $dag->name, $self->roi, $self->label));
        $self->__dag($dag);
    }
    return $self->__dag;
}

sub params_for_command {
    my $self = shift;
    return {
        label_fields => {
            roi_name => $self->roi,
            category => $self->label,
        },
        snvs_input_vcf => $self->input_vcf('snvs'),
        snvs_plan_file => $self->plan_file('snvs'),
        snvs_translations_file => $self->translations_file,
        indels_input_vcf => $self->input_vcf('indels'),
        indels_plan_file => $self->plan_file('indels'),
        indels_translations_file => $self->translations_file,
        use_header_from => 'snvs',
    };
}

sub roi {
    my $self = shift;
    return $self->discovery->region_of_interest_set->name;
}

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
    for my $report_plan ($plan->report_plans) {
        push @file_names, $report_plan->object->file_name;
    }
    return @file_names;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    my $translations_file = $self->generate_translations_file;
    $self->translations_file($translations_file);
    return $self;
};

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
        sprintf('Discovery(%s)', $self->discovery->tumor_sample->name)  => $self->discovery->tumor_bam,
        sprintf('Normal(%s)', $self->discovery->normal_sample->name) => $self->discovery->normal_bam,
    );
    if ($self->followup) {
        $bams{sprintf('Followup(%s)', $self->followup->tumor_sample->name)} = $self->followup->tumor_bam,
    }
    return %bams;
}

sub get_translations {
    my $self = shift;
    my %translations = %{$self->common_translations};
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

sub get_library_names {
    my $self = shift;

    my @instrument_data = $self->discovery->instrument_data;
    if (defined $self->followup) {
        push @instrument_data, $self->followup->instrument_data;
    }
    my @libraries = uniq map {$_->library} @instrument_data;
    return [map {$_->name} @libraries];
}

sub generate_translations_file {
    my $self = shift;

    my $translations = $self->get_translations;

    $translations->{aligned_bam_result_id} = $self->get_aligned_bams;
    $translations->{library_names} = $self->get_library_names;

    my %feature_list_ids;
    my $on_target_feature_list = Genome::FeatureList->get(name => $self->discovery->region_of_interest_set->name);
    $feature_list_ids{ON_TARGET} = $on_target_feature_list->id;
    my $segdup_feature_list = $self->discovery->get_feature_list("segmental_duplications");
    $feature_list_ids{SEG_DUP} = $segdup_feature_list->id;
    # TODO: There has to be a better way...
    $feature_list_ids{AML_RMG} = '0e4973c600244c3f804d54bee6f81145';
    $translations->{feature_list_ids} = \%feature_list_ids;
    $translations->{homopolymer_list_id} = '7f05e8fad6b6465a9f5bd6155dc88135';

    $translations->{reference_fasta} = $self->reference_sequence_build->full_consensus_path("fa");

    $translations->{dbsnp_vcf} = $self->discovery->previously_discovered_variations_build->snvs_vcf;
    $translations->{nhlbi_vcf} = _get_nhlbi_vcf(); 

    my $temp_file = Genome::Sys->create_temp_file_path;
    YAML::DumpFile($temp_file, $translations);

    return $temp_file;
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

