package Genome::VariantReporting::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);
use Set::Scalar;
use List::MoreUtils qw(uniq);
use File::Slurp qw();

my $DOCM = {
    snvs_build => "847b3cacad1249b8b7e46f89e02d96da",
    indels_build => "1040bf09070c4176a5256fa8a075378f",
};

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
        output_directory => {
            is => 'Path',
            is_output => 1,
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
    has_transient_optional => [
        _workflow => {
            is => 'Genome::WorkflowBuilder::DAG',
        },
        _workflow_inputs => {
            is => 'Hash',
        },
    ],
};

sub process_class {
    return "Genome::VariantReporting::Process::Trio";
}

sub execute {
    my $self = shift;

    my $process = $self->process_class->create();
    $self->initialize_dag;
    $self->add_summary_stats_to_dag;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        $self->add_reports_to_workflow($model_pair);
    }
    $self->add_merge_discovery_and_followup_reports_to_workflow;
    $self->add_igv_xml_to_workflow(\@model_pairs);
    $process->run(
        workflow_xml => $self->_workflow->get_xml
        workflow_inputs => $self->_workflow_inputs,
    );
    return 1;
}

sub add_igv_xml_to_workflow {
    my $self = shift;
    my $model_pairs = shift;

    my %bams = map { $_->get_sample_and_bam_map } @$model_pairs;
    my @reference_sequence_builds = uniq map { $_->reference_sequence_build } @$model_pairs;
    unless (scalar(@reference_sequence_builds) == 1) {
        die $self->error_message("Found more than one reference sequence build:" . Data::Dumper::Dumper(@reference_sequence_builds));
    }

    #TODO: find the right way to get these
    my @roi_directories = map {basename $_} glob(File::Spec->join($self->output_directory, "discovery", "*"));
    for my $roi_directory (@roi_directories) {

        my $params = {
            bam_hash_json => $_JSON_CODEC->canonical->encode(\%bams),
            genome_name => $self->tumor_sample->name,
            reference_name => $reference_sequence_builds[0]->name,
        };

        my $igv_op = Genome::WorkflowBuilder::Command->create(
            name => 'Create IGV session',
            command => 'Genome::VariantReporting::Command::Wrappers::CreateIGVSession'
        );

        $self->_add_operation_to_workflow(
            $igv_op,
            $params,
            "igv-$roi_name",
        );

        for my $bed_report_type (qw(discovery followup germline docm)) {
            $self->link_bed_report_to_igv($bed_report_type, $roi);
        }
    }
}

sub link_bed_report_to_igv {
    my $self = shift;
    #TODO
}

sub add_merge_discovery_and_followup_reports_to_workflow {
    my $self = shift;

    #TODO - find a better way of getting these
    my @roi_directories = map {basename $_} glob(File::Spec->join($self->output_directory, "discovery", "*"));
    for my $roi_directory (@roi_directories) {
        for my $base ($self->_trio_report_file_names) {
            my $discovery_report = File::Spec->join($self->output_directory, "discovery", $roi_directory, $base);
            my $additional_report = File::Spec->join($self->output_directory, "followup", $roi_directory, $base);
            my $merge_op = Genome::WorkflowBuilder::Command->new(
                name => "Merge discovery and followup reports ($report_name)",
                command => "Genome::VariantReporting::Framework::MergeReports",
            );
            my $params = {
                reports => [$discovery_report, $additional_report],
                sort_columns => [qw(chromosome_name start stop reference variant)],
                contains_header => 1,
                entry_sources =>  {$discovery_report => $self->tumor_sample->name, $additional_report => $self->followup_sample->name},
            };
        }
    }
    return 1;
}

sub _trio_report_file_names {
    return qw(trio_full_report.tsv trio_simple_report.tsv trio_acmg_report.tsv);
}

sub add_final_converge {
    my $self = shift;
    my $converge = Genome::WorkflowBuilder::Converge->create(
        name => "Final Converge",
        output_properties => ['reports'],
    );
    $self->_workflow->add_operation($converge);
    $self->_workflow->connect_output(
        output_property => "reports",
        source => $converge,
        source_property => "reports",
    );
}

sub get_model_pairs {
    my $self = shift;
    my $factory = Genome::VariantReporting::Command::Wrappers::ModelPairFactory->create(
        models => [$self->models],
        discovery_sample => $self->tumor_sample,
        followup_sample => $self->followup_sample,
        normal_sample => $self->normal_sample,
        output_dir => $self->output_directory,
        other_input_vcf_pairs => {docm => $self->vcf_files_from_imported_variation_builds($DOCM)},
    );
    return $factory->get_model_pairs;
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
    my ($self, $model_pair) = @_;

    $self->_add_operation_to_workflow(
        $model_pair->dag,
        $model_pair->params_for_dag,
        $model_pair->label
    );
}

sub initialize_dag {
    my $self = shift;
    my $dag = Genome::WorkflowBuilder::DAG->create(
        name => 'Trio Report',
        log_dir => File::Spec->join($self->output_directory, "logs_main"),
    );
    Genome::Sys->create_directory(File::Spec->join($self->output_directory, "logs_main"));
    $self->_workflow($dag);
    $self->add_final_converge;
}

sub _add_operation_to_workflow {
    my ($self, $operation, $params, $counter) = @_;
    $self->_workflow->add_operation($operation);
    for my $param (keys %$params) {
        my $full_param_name = $param;
        if (defined $counter) {
            $full_param_name = join("_", $param, $counter);
        }
        $self->_workflow->connect_input(
            input_property => $full_param_name,
            destination => $operation,
            destination_property => $param,
        );
        my $inputs = $self->_workflow_inputs;
        $inputs->{$full_param_name} = $params->{$param};
        $self->_workflow_inputs($inputs);
    }
}

sub add_summary_stats_to_dag {
    my $self = shift;
    if ($self->coverage_models) {
        my $alignment_stats_op = Genome::WorkflowBuilder::Command->create(
            name => 'Alignment stats',
            command => 'Genome::Model::SomaticValidation::Command::RunAlignmentStatsSummary'
        );
        my $coverage_stats_op = Genome::WorkflowBuilder::Command->create(
            name => 'Coverage stats',
            command => 'Genome::Model::SomaticValidation::Command::RunCoverageStatsSummary'
        );
        my $stats_params = {
            models => [$self->coverage_models],
        };
        $self->_add_operation_to_workflow(
            $alignment_stats_op,
            $stats_params,
        );
        $self->_add_operation_to_workflow(
            $coverage_stats_op,
            $stats_params,
        );
    }
}

1;

