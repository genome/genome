package Genome::VariantReporting::Command::Wrappers::Trio;

use strict;
use warnings;
use Genome;
use File::Basename qw(basename);
use Set::Scalar;
use List::MoreUtils qw(uniq);
use File::Slurp qw();

my $DOCM = {
    snvs_build => "a06152c107884ba78b90c5f09be17163",
    indels_build => "45841689d047419ea26df89128d6f121",
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

sub execute {
    my $self = shift;
    $self->initialize_dag;
    $self->run_summary_stats;
    my @model_pairs = $self->get_model_pairs;
    for my $model_pair (@model_pairs) {
        $self->add_reports_to_workflow($model_pair);
    }
    File::Slurp::write_file(File::Spec->join($self->output_directory, "workflow.xml"), $self->_workflow->get_xml);
    $self->_workflow->execute(%{$self->_workflow_inputs});
    $self->combine_discovery_and_followup_reports;
    $self->create_igv_xml(\@model_pairs);
    return 1;
}

sub create_igv_xml {
    my $self = shift;
    my $model_pairs = shift;

    my %bams = map { $_->get_sample_and_bam_map } @$model_pairs;
    my @reference_sequence_builds = uniq map { $_->reference_sequence_build } @$model_pairs;
    unless (scalar(@reference_sequence_builds) == 1) {
        die $self->error_message("Found more than one reference sequence build:" . Data::Dumper::Dumper(@reference_sequence_builds));
    }

    my @roi_directories = map {basename $_} glob(File::Spec->join($self->output_directory, "discovery", "*"));
    for my $roi_directory (@roi_directories) {
        my $discovery_bed = File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $self->output_directory, 'discovery', $roi_directory, 'trio_report.bed');
        my $additional_bed = File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $self->output_directory, 'followup', $roi_directory, 'trio_report.bed');
        my $germline_bed = File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $self->output_directory, 'germline', $roi_directory, 'germline_report.bed');
        my $docm_bed = File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $self->output_directory, 'docm', $roi_directory, 'cle_docm_report.bed');

        my $reference_sequence_name_cmd = Genome::Model::Tools::Analysis::ResolveIgvReferenceName->execute(
            reference_name => $reference_sequence_builds[0]->name,
        );

        #create the xml file for review
        my $dumpXML = Genome::Model::Tools::Analysis::DumpIgvXmlMulti->create(
            bams             => join(',', map {File::Spec->join($ENV{GENOME_SYS_SERVICES_FILES_URL}, $_)} values %bams),
            labels           => join(',', keys %bams),
            output_file      => File::Spec->join($self->output_directory, "$roi_directory.igv.xml"),
            genome_name      => $self->tumor_sample->name,
            review_bed_files => [$discovery_bed, $additional_bed, $germline_bed, $docm_bed],
            reference_name   => $reference_sequence_name_cmd->igv_reference_name,
        );
        unless ($dumpXML->execute) {
            confess $self->error_message("Failed to create IGV xml file");
        }
    }
}

sub combine_discovery_and_followup_reports {
    my $self = shift;

    my @roi_directories = map {basename $_} glob(File::Spec->join($self->output_directory, "discovery", "*"));
    for my $roi_directory (@roi_directories) {
        for my $base ($self->_trio_report_file_names) {
            my $discovery_report = File::Spec->join($self->output_directory, "discovery", $roi_directory, $base);
            my $additional_report = File::Spec->join($self->output_directory, "followup", $roi_directory, $base);
            Genome::VariantReporting::Command::CombineReports->execute(
                reports => [$discovery_report, $additional_report],
                sort_columns => [qw(chromosome_name start stop reference variant)],
                contains_header => 1,
                output_file => File::Spec->join($self->output_directory, "$base-$roi_directory"),
                entry_sources =>  {$discovery_report => $self->tumor_sample->name, $additional_report => $self->followup_sample->name},
            );
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
        $indels_build->snvs_vcf,
    ];
}

sub add_reports_to_workflow {
    my ($self, $model_pair) = @_;

    my %report_operations;
    for my $variant_type(qw(snvs indels)) {
        my %params = (
            input_vcf => $model_pair->input_vcf($variant_type),
            variant_type => $variant_type,
            output_directory => $model_pair->reports_directory($variant_type),
            plan_file => $model_pair->plan_file($variant_type),
            translations_file => $model_pair->translations_file,
            log_directory => $model_pair->logs_directory($variant_type),
        );
        $report_operations{$variant_type} = $self->add_report_to_workflow(\%params);
    }
    my $snv_reports = Set::Scalar->new(grep {!($_ =~ /vcf$/)} $model_pair->report_names("snvs"));
    my $indel_reports = Set::Scalar->new(grep {!($_ =~ /vcf$/)} $model_pair->report_names("indels"));
    my $both_reports = $snv_reports->intersection($indel_reports);
    for my $base ($both_reports->members) {
        my %combine_params = (
            output_file => File::Spec->join($model_pair->output_dir, $base),
            separator => "\t",
        );
        if ($base =~ m/bed$/) {
            $combine_params{sort_columns} = [qw(1 2 3)];
            $combine_params{contains_header} = 0,
        }
        else {
            $combine_params{sort_columns} = [qw(chromosome_name start stop reference variant)];
            $combine_params{contains_header} = 1,
            $combine_params{split_indicators} = [qw(per_library)],
        }
        $self->add_combine_to_workflow(\%combine_params, \%report_operations, $base);
    }
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

sub add_combine_to_workflow {
    my ($self, $params, $reports_to_combine, $file_name) = @_;
    #snv report and indel report => combine_snv_indel_reports helper => combine_reports => output_connector
    my $counter = $self->_workflow_counter;
    my $combine_command = Genome::WorkflowBuilder::Command->create(
        name => join(" ", "Combine snv and indel", $counter),
        command => "Genome::VariantReporting::Command::CombineReports",
    );
    $self->_add_operation_to_workflow($combine_command, $params, $counter);
    my $generate_file_names = Genome::WorkflowBuilder::Command->create(
        name => join(" ", "Generate filenames", $counter),
        command => "Genome::VariantReporting::Command::CombineSnvIndelReports",
    );
    $self->_workflow->add_operation($generate_file_names);
    my $converge = Genome::WorkflowBuilder::Converge->create(
        output_properties => ['output_dirs'],
        name => join(" ", "Converge", $counter),
    );
    $self->_workflow->add_operation($converge);
    while (my ($variant_type, $report) = each %$reports_to_combine) {
        $self->_workflow->create_link(
            source => $report, source_property => "output_directory",
            destination => $converge, destination_property => $variant_type."_output_dir",
        );
    }
    my $full_param_name = join("_", "file_name", $counter);
    $self->_workflow->connect_input(
            input_property => $full_param_name,
            destination => $generate_file_names,
            destination_property => "file_name",
    );
    $self->_workflow_inputs->{$full_param_name} = $file_name;
    $self->_workflow->create_link(
        source => $converge, source_property => "output_dirs",
        destination => $generate_file_names, destination_property => "input_directories",
    );
    $self->_workflow->create_link(
        source => $generate_file_names, source_property => "reports",
        destination => $combine_command, destination_property => "reports"
    );
    my $final_converge = $self->_workflow->operation_named("Final Converge");
    $self->_workflow->create_link(
        source => $combine_command, source_property => "output_file",
        destination => $final_converge, destination_property => $combine_command->name,
    );
}

sub add_report_to_workflow {
    my ($self, $params) = @_;
    my $counter = $self->_workflow_counter;
    my $report_creator = Genome::VariantReporting::Command::CreateReport->create(%$params);
    my $report_dag = $report_creator->dag;
    $report_dag->name(join(" ", $report_dag->name, $counter));
    $self->_add_operation_to_workflow($report_dag, {$report_creator->params_for_execute}, $counter);
    return $report_dag;
}

sub _add_operation_to_workflow {
    my ($self, $operation, $params, $counter) = @_;
    $self->_workflow->add_operation($operation);
    for my $param (keys %$params) {
        my $full_param_name = join("_", $param, $counter);
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

sub run_summary_stats {
    my $self = shift;
    if ($self->coverage_models) {
        Genome::Model::SomaticValidation::Command::AlignmentStatsSummary->execute(
            output_tsv_file => File::Spec->join($self->output_directory, "alignment_summary.tsv"),
            models => [$self->coverage_models],
        );
        Genome::Model::SomaticValidation::Command::CoverageStatsSummary->execute(
            output_tsv_file => File::Spec->join($self->output_directory, "coverage_summary.tsv"),
            models => [$self->coverage_models],
        );
    }
}

{
    my $counter = 0;
    sub _workflow_counter {
        return $counter++;
    }
}

1;

