package Genome::Model::Tools::EpitopePrediction::NewPipeline;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::EpitopePrediction::NewPipeline {
    is => 'Command::V2',
    has_input => [
        output_directory => {
            is => 'Text',
            doc => 'the directory where you want results stored',
        },
        build => {
            is => "Genome::Model::Build::SomaticVariation",
            doc => 'The somatic variation build whose variants are to be used for analysis',
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

    Genome::Sys->create_directory($self->output_directory);
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

    unless (defined($self->{dag})) {
        my $dag = Genome::WorkflowBuilder::DAG->create(
            name => 'Epitope Binding Predicition',
        );

        $self->add_reports_to_workflow($dag);
        $self->{dag} = $dag;
    }

    return $self->{dag};
}

sub add_reports_to_workflow {
    my ($self, $dag) = @_;

    for my $variant_type (qw(snvs indels)) {
        my %report_params = (
            build => $self->build,
            variant_type => $variant_type,
        );
        my $wrapper = Genome::VariantReporting::Command::Wrappers::EpitopeBindingPredictionFasta->create(%report_params);
        next unless $wrapper->has_valid_variant_type_for_build;

        my $report_dag = $wrapper->report_workflow;
        $dag->add_operation($report_dag);
        $dag->connect_input(
            input_property => 'process_id',
            destination => $report_dag,
            destination_property => 'process_id',
        );
        $wrapper->delete();

        my $output_directory = File::Spec->join($self->output_directory, $variant_type);
        Genome::Sys->create_directory($output_directory);
        my %netmhc_params = (
            output_directory => $output_directory,
            epitope_length => $self->epitope_length,
            netmhc_version => $self->netmhc_version,
            output_filter => $self->output_filter,
            sample_name =>  $self->build->subject_name,
            alleles => [$self->alleles],
            variant_type => $variant_type,
        );
        my $netmhc = Genome::Model::Tools::EpitopePrediction::NetmhcPipeline->create(%netmhc_params);
        my $netmhc_dag = $netmhc->_construct_workflow;
        my %netmhc_inputs = %{$netmhc->_get_workflow_inputs};
        delete $netmhc_inputs{input_fasta_file};
        $netmhc_dag->declare_constant(%netmhc_inputs);
        $dag->add_operation($netmhc_dag);
        $netmhc->delete();

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

1;
