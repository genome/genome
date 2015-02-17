package Genome::VariantReporting::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);
use Set::Scalar;
use List::MoreUtils qw(uniq);
use File::Slurp qw();
use JSON qw(to_json);

my $DOCM = {
    snvs_build => "847b3cacad1249b8b7e46f89e02d96da",
    indels_build => "1040bf09070c4176a5256fa8a075378f",
};

my $_JSON_CODEC = new JSON->allow_nonref;
my $FACTORY = Genome::VariantReporting::Framework::Factory->create();

class Genome::VariantReporting::Command::Wrappers::Trio {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            doc => "Models to run variant reports on",
        },
        coverage_models => {
            is => 'Genome::Model::SomaticValidation',
            is_many => 1,
            is_optional => 1,
            doc => "Models to run coverage reports on",
        },
        tumor_sample => {
            is => 'Genome::Sample',
            doc => 'Main tumor sample used for discovery',
        },
        followup_sample => {
            is => 'Genome::Sample',
            doc => 'Additional sample to report readcounts on at discovery variant positions',
        },
        normal_sample => {
            is => 'Genome::Sample',
            doc => 'Normal sample',
        },
    ],
};

sub process_class {
    return "Genome::VariantReporting::Process::Trio";
}

sub execute {
    my $self = shift;

    my $p = $self->process_class->create(
        builds => [$self->builds],
        coverage_builds => [$self->coverage_builds],
        tumor_sample => $self->tumor_sample,
        followup_sample => $self->followup_sample,
        normal_sample => $self->normal_sample,
    );

    $self->status_message("Constructing workflow from inputs. (this may take a while...)");
    $p->run(workflow_xml => $self->dag->get_xml,
        workflow_inputs => $self->workflow_inputs($p->id),
    );

    return $p;
}

sub builds {
    my $self = shift;
    return map {$_->last_succeeded_build} $self->models;
}

sub coverage_builds {
    my $self = shift;
    return map {$_->last_succeeded_build} $self->coverage_models;
}

sub workflow_inputs {
    my $self = shift;
    my $process_id = shift;
    return {
        process_id => $process_id,
        %{$self->dag->constant_values},
    };
}

sub dag {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'Trio Report',
    );

    $self->add_summary_stats_to_dag($dag);
    my ($model_pairs, $models_for_roi) = $self->get_model_pairs_and_models_for_roi;

    for my $model_pair (@{$model_pairs}) {
        $self->add_reports_to_workflow($dag, $model_pair);
    }

    $self->add_merge_discovery_and_followup_reports_to_workflow(
        $dag, keys %{$models_for_roi});
    $self->add_igv_xml_to_workflow($dag, $model_pairs, keys %{$models_for_roi});

    return $dag;
}

sub add_igv_xml_to_workflow {
    my $self = shift;
    my $dag = shift;
    my $model_pairs = shift;
    my @roi_names = @_;

    my %bams = map { $_->get_sample_and_bam_map } @$model_pairs;
    my @reference_sequence_builds = uniq map { $_->reference_sequence_build } @$model_pairs;
    unless (scalar(@reference_sequence_builds) == 1) {
        die $self->error_message("Found more than one reference sequence build:" .
            Data::Dumper::Dumper(@reference_sequence_builds));
    }

    for my $roi_name (@roi_names) {
        my $converge = Genome::WorkflowBuilder::Converge->create(
            output_properties => ['report_results'],
            name => "Converge for igv ($roi_name)",
        );
        $dag->add_operation($converge);

        for my $category (qw(discovery followup germline), keys %{$self->other_input_vcf_pairs}) {
            my $sub_dag = $dag->operation_named(sub_dag_name($roi_name, $category));
            my $output_name = "merged_result (bed)";

            $dag->create_link(
                source => $sub_dag,
                source_property => $output_name,
                destination => $converge,
                destination_property => sprintf('%s_bed_report_result', $category),
            );

        }
        my $igv_op = Genome::WorkflowBuilder::Command->create(
            name => sprintf('Create IGV session (%s)', $roi_name),
            command => 'Genome::VariantReporting::Command::Wrappers::CreateIgvSession'
        );
        $dag->create_link(
            source => $converge,
            source_property => 'report_results',
            destination => $igv_op,
            destination_property => 'merged_bed_reports',
        );
        $igv_op->declare_constant(
            bam_hash_json => $_JSON_CODEC->canonical->encode(\%bams),
            genome_name => $self->tumor_sample->name,
            reference_name => $reference_sequence_builds[0]->name,
            label => 'igv_session:' . to_json({
                'roi_name' => $roi_name
            }, {canonical => 1}),
        );
        $dag->connect_input(
            input_property => 'process_id',
            destination => $igv_op,
            destination_property => 'process_id',
        );

        $dag->connect_output(
            output_property => sprintf("igv_session.%s", $roi_name),
            source => $igv_op,
            source_property => 'output_result',
        );
        $dag->add_operation($igv_op);
    }
    return;
}

sub link_bed_report_to_igv {
    my $self = shift;
    #TODO
}

sub sub_dag_name {
    my ($roi_name, $category) = @_;

    # TODO do something less smart.
    return sprintf('Create Snvs, Indels, and Merged Reports (%s-%s)',
        $roi_name, $category);
}

sub add_merge_discovery_and_followup_reports_to_workflow {
    my $self = shift;
    my $dag = shift;
    my @roi_names = @_;

    for my $roi_name (@roi_names) {
        my $discovery_dag = $dag->operation_named(
            sub_dag_name($roi_name, 'discovery'));
        for my $output_name ($discovery_dag->output_properties) {
            if ($output_name =~ m/merged_result \((.*)\)/) {
                my $report_name = $1;
                my $report_class = $FACTORY->get_class('reports', $report_name);
                next unless $report_class->can_be_merged;

                my $merge_op = Genome::WorkflowBuilder::Command->create(
                    name => sprintf('Merge Discovery and Followup Reports (%s)', $report_name),
                    command => 'Genome::VariantReporting::Framework::MergeReports',
                );

                $dag->create_link(
                    source => $discovery_dag,
                    source_property => $output_name,
                    destination => $merge_op,
                    destination_property => 'base_report',
                );

                my $followup_dag = $dag->operation_named(sub_dag_name(
                        $roi_name, 'followup'));
                $dag->create_link(
                    source => $followup_dag,
                    source_property => $output_name,
                    destination => $merge_op,
                    destination_property => 'supplemental_report',
                );

                $merge_op->declare_constant(
                    label => 'report:' . to_json({
                            roi_name => $roi_name,
                            category => 'discovery_and_followup',
                            variant_type => 'merged',
                            report_name => $report_name,
                        }, {canonical=>1}),
                    %{$report_class->merge_parameters},
                );
                # this has to be done AFTER the constants are declared.
                $dag->add_operation($merge_op);

                $dag->connect_input(
                    input_property => 'process_id',
                    destination => $merge_op,
                    destination_property => 'process_id',
                );

                $dag->connect_output(
                    output_property => sprintf('%s-discovery_and_followup.merged.%s',
                        $roi_name, $report_name),
                    source => $merge_op,
                    source_property => 'output_result',
                );
            }
        }
    }
    return;
}

sub _trio_report_file_names {
    return qw(trio_full_report.tsv trio_simple_report.tsv trio_acmg_report.tsv);
}

sub other_input_vcf_pairs {
    my $self = shift;

    return {docm => $self->vcf_files_from_imported_variation_builds($DOCM)};
}

sub get_model_pairs_and_models_for_roi {
    my $self = shift;
    my $factory = Genome::VariantReporting::Command::Wrappers::ModelPairFactory->create(
        models => [$self->models],
        discovery_sample => $self->tumor_sample,
        followup_sample => $self->followup_sample,
        normal_sample => $self->normal_sample,
        other_input_vcf_pairs => $self->other_input_vcf_pairs,
    );
    return $factory->get_model_pairs, $factory->get_models_for_roi;
}

sub vcf_files_from_imported_variation_builds {
    my ($self, $builds) = @_;
    my $snvs_build = Genome::Model::Build->get($builds->{snvs_build});
    my $indels_build = Genome::Model::Build->get($builds->{indels_build});
    return [
        $snvs_build->snvs_vcf,
        $indels_build->indels_vcf,
    ];
}

sub add_reports_to_workflow {
    my ($self, $dag, $model_pair) = @_;

    my $model_pair_dag = $model_pair->dag;
    $dag->add_operation($model_pair_dag);
    $dag->connect_input(
        input_property => 'process_id',
        destination => $model_pair_dag,
        destination_property => 'process_id',
    );

    return;
}

sub add_summary_stats_to_dag {
    my $self = shift;
    my $dag = shift;

    if ($self->coverage_models) {
        my %commands = (
                'Alignment stats' => 'Genome::Model::SomaticValidation::Command::RunAlignmentStatsSummary',
                'Coverage stats' => 'Genome::Model::SomaticValidation::Command::RunCoverageStatsSummary',
        );
        while (my($name, $command) = each %commands) {
            my $op = Genome::WorkflowBuilder::Command->create(
                name => $name,
                command => $command,
            );
            $op->declare_constant(
                models => [$self->coverage_models],
            );
            $dag->connect_input(
                input_property => 'process_id',
                destination => $op,
                destination_property => 'process_id',
            );
            $dag->connect_output(
                output_property => sprintf("%s result", $name),
                source => $op,
                source_property => 'output_result',
            );
            $dag->add_operation($op);
        }
    }
    return;
}

1;

