package Genome::Model::Tools::EpitopePrediction::NetmhcPipeline;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::EpitopePrediction::NetmhcPipeline {
    is => 'Command::V2',
    doc => 'Run the netmhc portion of the epitope binding prediction pipeline',
    has_input => [
        output_directory => {
            is => 'Text',
            doc => 'the directory where you want results stored',
        },
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_optional => 1,
            doc => 'The somatic variation build to use for analysis',
        },
        input_fasta_file => {
            is => 'Text',
            doc => 'The fasta file with variant sequences for epitope binding prediction',
        },
        alleles => {
            is => 'Text',
            doc => 'A list of allele names to be used for epitope prediction with NetMHC',
            is_many => 1,
        },
        epitope_length => {
            is => 'Integer',
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
        sample_name => {
            is => 'Text',
            doc => 'The sample name of the file being processed',
            is_optional => 1,
        },
        variant_type => {
            is => 'Text',
            doc => 'The variant type being processed',
            is_optional => 1,
        }
    ],
};

sub command_class_prefix {
    return "Genome::Model::Tools::EpitopePrediction";
}

sub execute {
    my $self = shift;

    $self->debug_message("Validating Inputs...");
    $self->_validate_inputs();

    $self->debug_message("Constructing Workflow...");
    my $workflow = $self->_construct_workflow();

    $self->debug_message("Getting Workflow Inputs...");
    my $inputs = $self->_get_workflow_inputs();

    $self->debug_message("Running Workflow...");
    my $result = $workflow->execute(inputs => $inputs);

    unless($result){
        $self->fatal_message("Workflow did not return correctly.");
    }

    return 1;
}

sub _construct_workflow {
    my ($self) = @_;

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => $self->_workflow_name,
        log_dir => $self->output_directory,
    );

    my $filter_sequences_command = $self->_attach_filter_sequences_command($workflow);
    my $generate_fasta_key_command = $self->_attach_generate_fasta_key_command($workflow);

    $workflow->create_link(
        source => $filter_sequences_command,
        source_property => 'output_file',
        destination => $generate_fasta_key_command,
        destination_property => 'input_file',
    );

    my $netmhc_workflow = $self->create_netmhc_workflow;
    $workflow->add_operation($netmhc_workflow);
    for my $property (qw/allele epitope_length netmhc_version sample_name output_directory output_filter/) {
        $workflow->connect_input(
            input_property => $property,
            destination => $netmhc_workflow,
            destination_property => $property,
        );
    }

    $workflow->create_link(
        source => $filter_sequences_command,
        source_property => 'output_file',
        destination => $netmhc_workflow,
        destination_property => 'fasta_file',
    );
    $workflow->create_link(
        source => $generate_fasta_key_command,
        source_property => 'output_file',
        destination => $netmhc_workflow,
        destination_property => 'key_file',
    );

    $workflow->connect_output(
        output_property => "output_file",
        source => $netmhc_workflow,
        source_property => 'output_file',
    );

    return $workflow;
}

sub _workflow_name {
    my $self = shift;

    if (defined($self->variant_type)) {
        return sprintf('Epitope Prediction Workflow (%s)', $self->variant_type);
    }
    else {
        return 'Epitope Prediction Workflow';
    }
}

sub create_netmhc_workflow {
    my $self = shift;

    my $netmhc_workflow = Genome::WorkflowBuilder::DAG->create(
        name => $self->_netmhc_workflow_name,
        log_dir => $self->output_directory,
        parallel_by => 'allele',
    );
    my $run_netmhc_command = $self->_attach_run_netmhc_command($netmhc_workflow);
    my $parse_netmhc_command = $self->_attach_parse_netmhc_command($netmhc_workflow);

    $netmhc_workflow->create_link(
        source => $run_netmhc_command,
        source_property => 'output_file',
        destination => $parse_netmhc_command,
        destination_property => 'netmhc_file',
    );

    $netmhc_workflow->connect_output(
        output_property => "output_file",
        source => $parse_netmhc_command,
        source_property => 'parsed_file',
    );

    return $netmhc_workflow;
}

sub _netmhc_workflow_name {
    my $self = shift;

    if (defined($self->variant_type)) {
        return sprintf('NetMHC Workflow (%s)', $self->variant_type);
    }
    else {
        return 'NetMHC Workflow';
    }

}

sub _attach_filter_sequences_command {
    my $self = shift;
    my $workflow = shift;

    my $filter_sequences_command = Genome::WorkflowBuilder::Command->create(
        name => 'FilterSequencesCommand',
        command => $self->filter_sequences_command_name,
    );
    $workflow->add_operation($filter_sequences_command);
    $self->_add_common_inputs($workflow, $filter_sequences_command);
    $workflow->connect_input(
        input_property => 'input_fasta_file',
        destination => $filter_sequences_command,
        destination_property => 'input_file',
    );
    return $filter_sequences_command;
}

sub _attach_generate_fasta_key_command {
    my $self = shift;
    my $workflow = shift;

    my $generate_fasta_key_command = Genome::WorkflowBuilder::Command->create(
        name => 'GenerateFastaKeyCommand',
        command => $self->generate_fasta_key_command_name,
    );
    $workflow->add_operation($generate_fasta_key_command);
    $self->_add_common_inputs($workflow, $generate_fasta_key_command);
    return $generate_fasta_key_command;
}

sub _attach_run_netmhc_command {
    my $self = shift;
    my $workflow = shift;

    my $run_netmhc_command = Genome::WorkflowBuilder::Command->create(
        name => "RunNetMHCCommand",
        command => $self->run_netmhc_command_name,
    );
    $workflow->add_operation($run_netmhc_command);
    for my $property (qw/allele epitope_length netmhc_version sample_name fasta_file output_directory/) {
        $workflow->connect_input(
            input_property => $property,
            destination => $run_netmhc_command,
            destination_property => $property,
        );
    }

    return $run_netmhc_command;
}

sub _attach_parse_netmhc_command{
    my $self = shift;
    my $workflow = shift;

    my $parse_netmhc_command = Genome::WorkflowBuilder::Command->create(
        name => "ParseNetMHCCommand",
        command => $self->parse_netmhc_command_name,
    );
    $workflow->add_operation($parse_netmhc_command);
    for my $property (qw/output_filter netmhc_version key_file output_directory/) {
        $workflow->connect_input(
            input_property => $property,
            destination => $parse_netmhc_command,
            destination_property => $property,
        );
    }
    return $parse_netmhc_command;
}

sub get_wildtype_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::GetWildtype";
}

sub generate_variant_sequences_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::GenerateVariantSequences";
}

sub filter_sequences_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::FilterSequences";
}

sub generate_fasta_key_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::GenerateFastaKey";
}

sub run_netmhc_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::RunNetmhc";
}

sub parse_netmhc_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::ParseNetmhcOutput";
}

sub _add_common_inputs {
    my $self = shift;
    my $workflow = shift;
    my $command = shift;

    my @common_inputs = qw(
        output_directory
    );

    for my $prop_name (@common_inputs) {
        $workflow->connect_input(
            input_property => $prop_name,
            destination => $command,
            destination_property => $prop_name,
        );
    }
}

sub _validate_inputs {
    my $self = shift;

    if (defined($self->somatic_variation_build)) {
        if (defined($self->sample_name)) {
            $self->status_message("Custom sample name provided. Using custom sample name %s instead of somatic variation build sample name", $self->sample_name);
        }
        else {
            my $sample_name = $self->somatic_variation_build->subject_name;
            $self->status_message("Somatic variation build given. Setting sample name to $sample_name");
            $self->sample_name($sample_name);
        }
    }
    else {
        unless (defined($self->sample_name)) {
            $self->fatal_message("Sample name must be defined if no somatic variation build is given")
        }
    }

    unless (-s $self->input_fasta_file) {
        $self->fatal_message("Input fasta file %s does not exist or has no size", $self->input_fasta_file);
    }

    unless (Genome::Sys->create_directory($self->output_directory)) {
        $self->fatal_message("Could not create directory (%s)", $self->output_directory);
    }

    for my $allele ($self->alleles) {
        unless (Genome::Model::Tools::EpitopePrediction::RunNetmhc->is_valid_allele_for_netmhc_version($allele, $self->netmhc_version)) {
            $self->fatal_message("Allele %s not valid for NetMHC version %s", $allele, $self->netmhc_version);
        }
    }

    return 1;
}

sub _get_workflow_inputs {
    my $self = shift;

    my %inputs = (
        input_fasta_file => $self->input_fasta_file,
        output_directory => $self->output_directory,
        epitope_length => $self->epitope_length,
        netmhc_version => $self->netmhc_version,
        output_filter => $self->output_filter,
        sample_name => $self->sample_name,
        allele => [$self->alleles],
    );

    return \%inputs;
}

sub final_output_file {
    my $self = shift;
    my $allele = shift;

    my $file_name = join ('.', $self->sample_name, $allele, $self->epitope_length, 'netmhc', 'parsed', $self->output_filter);
    return File::Spec->join($self->output_directory, $file_name);
}

1;
