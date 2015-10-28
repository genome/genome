package Genome::Model::Tools::EpitopePrediction::NewPipeline;

use strict;
use warnings;
use Genome;
use File::Basename qw(dirname);
use YAML;

class Genome::Model::Tools::EpitopePrediction::NewPipeline {
    is => 'Command::V2',
    has_input => [
        output_directory => {
            is => 'Text',
            doc => 'the directory where you want results stored',
        },
        build => {
            is => "Genome::Model::Build::SomaticVariation",
        },
        alleles => {
            is => 'Text',
            doc => 'A list of allele names to be used for epitope prediction with NetMHC',
            is_many => 1,
        },
        epitope_length => {
            is => 'Text',
            doc => 'Length of subpeptides to predict with NetMHC',
        },
        netmhc_version => {
            is => 'Text',
            doc => 'The NetMHC version to use',
            valid_values => ['3.0','3.4'],
            default_value => '3.4',
        },
        output_filter => {
            is => 'Text',
            doc =>
                'Type of epitopes to report in the final output - select \'top\' to report the top epitopes in terms of fold changes,  \'all\' to report all predictions ',
            valid_values => ['top', 'all'],
        },
    ],
};

sub process_class {
    return "Genome::Model::Tools::EpitopePrediction::Process";
}

sub execute {
    my $self = shift;

    $self->validate_inputs;

    my $p = $self->process_class->create(
        $self->process_params,
    );

    $self->status_message("Constructing workflow from inputs. (this may take a while...)");
    $p->run(
        workflow_xml => $self->dag->get_xml,
        workflow_inputs => $self->workflow_inputs($p->id),
    );

    return $p;
}

sub validate_inputs {
    my $self = shift;

    if (-d $self->output_directory) {
        $self->fatal_message("Output directory %s already exists", $self->output_directory);
    }

    Genome::Sys::make_path($self->output_directory);
}

sub process_params {
    my $self = shift;
    return (
        build => $self->build,
    );
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
        name => 'Epitope Binding Predicition',
    );

    $self->add_reports_to_workflow($dag);

    return $dag;
}

sub add_reports_to_workflow {
    my ($self, $dag) = @_;

    for my $variant_type (qw(snvs indels)) {
        my $input_vcf = $self->input_vcf($variant_type);
        next unless -s $input_vcf;

        my %report_params = (
            input_vcf => $input_vcf,
            variant_type => $variant_type,
            plan_file => $self->plan_file($variant_type),
            translations_file => $self->translations_file,
        );
        my $report = Genome::VariantReporting::Command::CreateReport->create(%report_params);
        my $report_dag = $report->dag;
        $dag->add_operation($report_dag);
        $dag->connect_input(
            input_property => 'process_id',
            destination => $report_dag,
            destination_property => 'process_id',
        );

        my $output_directory = File::Spec->join($self->output_directory, $variant_type);
        Genome::Sys::make_path($output_directory);
        my %netmhc_params = (
            output_directory => $output_directory,
            epitope_length => $self->epitope_length,
            netmhc_version => $self->netmhc_version,
            output_filter => $self->output_filter,
            sample_name =>  $self->build->subject_name,
            alleles => $self->alleles,
            variant_type => $variant_type,
        );
        my $netmhc = Genome::Model::Tools::EpitopePrediction::NetmhcPipeline->create(%netmhc_params);
        my $netmhc_dag = $netmhc->_construct_workflow;
        my %netmhc_inputs = %{$netmhc->_get_workflow_inputs};
        delete $netmhc_inputs{input_fasta_file};
        $netmhc_dag->declare_constant(%netmhc_inputs);
        $dag->add_operation($netmhc_dag);

        for my $report_output_name ($report_dag->output_properties) {
            if ($report_output_name =~ m/report_path \((.*)\)/) {
                my $report_name = $1;
                $dag->create_link(
                    source => $report_dag,
                    source_property => $report_output_name,
                    destination => $netmhc_dag,
                    destination_property => 'input_fasta_file',
                );
                $dag->connect_output(
                    output_property => sprintf('%s_result (%s)', $variant_type, $report_name),
                    source => $netmhc_dag,
                    source_property => 'output_file',
                );
            }
        }
    }

    return;
}

sub translations_file {
    my $self = shift;
    my $translations_file = Genome::Sys->create_temp_file_path;
    my %translations;
    $translations{reference_fasta} = $self->build->reference_sequence_build->full_consensus_path("fa");
    YAML::DumpFile(File::Spec->join($translations_file), \%translations);
    return $translations_file;
}

sub plan_file {
    my ($self, $variant_type) = @_;
    my $base_dir = dirname(dirname(dirname(dirname(__FILE__))));
    return File::Spec->join($base_dir, 'VariantReporting', 'plan_files', "epitope_prediction_$variant_type.yaml");
}

sub input_vcf {
    my ($self, $variant_type) = @_;
    my $accessor = "get_detailed_${variant_type}_vcf";
    return $self->build->$accessor;
}

1;
