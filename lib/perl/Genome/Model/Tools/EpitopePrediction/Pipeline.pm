package Genome::Model::Tools::EpitopePrediction::Pipeline;

use strict;
use warnings;

use Genome;
use Workflow::Simple;

class Genome::Model::Tools::EpitopePrediction::Pipeline {
    is => 'Command::V2',
    doc => 'Run the epitope binding prediction pipeline',
    has => [
        output_directory => {
            is => 'Text',
            doc => 'the directory where you want results stored',
        },
        #TODO fill out docs
        somatic_variation_build => {
            is => 'Genome::Model::Build::SomaticVariation',
            is_optional => 1,
            doc => '',
        },
        input_tsv_file => {
            is => 'Text',
            is_optional => 1,
            doc => '',
        },
        anno_db => {
            is => 'Text',
            is_optional => 1,
            doc => 'The name of the annotation database to use for retrieving the wildtypes.  Example: NCBI-human.combined-annotation',
        },
        anno_db_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'The version of the annotation databaseto use for retrieving the wildtypes. Example: 54_36p_v2',
        },
        peptide_sequence_length => {
            is => 'Text',
            doc => 'The length of the peptide sequences to be used when generating variant sequences',
            valid_values => [17, 21, 31],
            default_value => 21,
        },
        allele => {
            is => 'Text',
            doc => 'Allele name to be used for epitope prediction with NetMHC',
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
        sample_name => {
            is => 'Text',
            doc => 'The sample name of the file being processed',
            is_optional => 1,
        },
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
    my $result = $workflow->execute(%$inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

sub _construct_workflow {
    my ($self) = @_;

    # my $xml = __FILE__ . '.xml';
    # my $workflow = Workflow::Operation->create_from_xml($xml);
    # $workflow->log_dir($self->output_directory);

    my $workflow = Genome::WorkflowBuilder::DAG->create(
        name => 'EpitopePredictionWorkflow',
        log_dir => $self->output_directory,
    );

    my $get_wildtype_command = $self->_attach_get_wildtype_command($workflow);
    my $generate_variant_sequences_command = $self->_attach_generate_variant_sequences_command($workflow);
    my $remove_star_sequences_command = $self->_attach_remove_star_sequences_command($workflow);
    my $generate_fasta_key_command = $self->_attach_generate_fasta_key_command($workflow);
    my $run_netmhc_command = $self->_attach_run_netmhc_command($workflow);
    my $parse_netmhc_command = $self->_attach_parse_netmhc_command($workflow);

    $workflow->create_link(
        source => $get_wildtype_command,
        source_property => 'output_tsv_file',
        destination => $generate_variant_sequences_command,
        destination_property => 'input_file',
    );
    $workflow->create_link(
        source => $generate_variant_sequences_command,
        source_property => 'output_file',
        destination => $remove_star_sequences_command,
        destination_property => 'input_file',
    );
    $workflow->create_link(
        source => $remove_star_sequences_command,
        source_property => 'output_file',
        destination => $generate_fasta_key_command,
        destination_property => 'input_file',
    );
    $workflow->create_link(
        source => $remove_star_sequences_command,
        source_property => 'output_file',
        destination => $run_netmhc_command,
        destination_property => 'fasta_file',
    );
    $workflow->create_link(
        source => $run_netmhc_command,
        source_property => 'output_file',
        destination => $parse_netmhc_command,
        destination_property => 'netmhc_file',
    );
    $workflow->create_link(
        source => $generate_fasta_key_command,
        source_property => 'output_file',
        destination => $parse_netmhc_command,
        destination_property => 'key_file',
    );

    $workflow->connect_output(
        output_property => 'output_file',
        source => $parse_netmhc_command,
        source_property => 'parsed_file',
    );

    return $workflow;
}

sub _attach_get_wildtype_command {
    my $self = shift;
    my $workflow = shift;

    my $get_wildtype_command = Genome::WorkflowBuilder::Command->create(
        name => 'GetWildTypeCommand',
        command => $self->get_wildtype_command_name,
    );
    $workflow->add_operation($get_wildtype_command);
    $self->_add_common_inputs($workflow, $get_wildtype_command);
    for my $property (qw/input_tsv_file anno_db anno_db_version/) {
        $workflow->connect_input(
            input_property => $property,
            destination => $get_wildtype_command,
            destination_property => $property,
        );
    }
    return $get_wildtype_command;
}

sub _attach_generate_variant_sequences_command {
    my $self = shift;
    my $workflow = shift;

    my $generate_variant_sequences_command = Genome::WorkflowBuilder::Command->create(
        name => 'GenerateVariantSequencesCommand',
        command => $self->generate_variant_sequences_command_name,
    );
    $workflow->add_operation($generate_variant_sequences_command);
    $self->_add_common_inputs($workflow, $generate_variant_sequences_command);
    for my $property (qw/peptide_sequence_length/) {
        $workflow->connect_input(
            input_property => $property,
            destination => $generate_variant_sequences_command,
            destination_property => $property,
        );
    }
    return $generate_variant_sequences_command;
}

sub _attach_remove_star_sequences_command {
    my $self = shift;
    my $workflow = shift;

    my $remove_star_sequences_command = Genome::WorkflowBuilder::Command->create(
        name => 'RemoveStarSequencesCommand',
        command => $self->remove_star_sequences_command_name,
    );
    $workflow->add_operation($remove_star_sequences_command);
    $self->_add_common_inputs($workflow, $remove_star_sequences_command);
    return $remove_star_sequences_command;
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
        name => 'RunNetMHCCommand',
        command => $self->run_netmhc_command_name,
    );
    $workflow->add_operation($run_netmhc_command);
    $self->_add_common_inputs($workflow, $run_netmhc_command);
    for my $property (qw/allele epitope_length netmhc_version sample_name/) {
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
        name => 'ParseNetMHCCommand',
        command => $self->parse_netmhc_command_name,
    );
    $workflow->add_operation($parse_netmhc_command);
    $self->_add_common_inputs($workflow, $parse_netmhc_command);
    for my $property (qw/output_filter netmhc_version/) {
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

sub remove_star_sequences_command_name {
    my $self = shift;

    return $self->command_class_prefix . "::RemoveStarSequences";
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

    if (!defined($self->somatic_variation_build) && !defined($self->input_tsv_file)) {
        die $self->error_message("Either somatic variation build or input tsv file needs to be provided");
    }

    if (defined($self->somatic_variation_build)) {
        if (defined($self->input_tsv_file)) {
            die $self->error_message("Custom tsv file cannot be used in combination with somatic variation build");
        }
        else {
            my $tsv_file = File::Spec->join(
                $self->somatic_variation_build->data_directory,
                'effects',
                'snvs.hq.tier1.v1.annotated.top.header'
            );
            $self->status_message("Somatic variation build given. Setting input_tsv_file to $tsv_file");
            $self->input_tsv_file($tsv_file);
        }

        if (defined($self->sample_name)) {
            die $self->error_message("Custom sample name cannot be used in combination with somatic variation build");
        }
        else {
            my $sample_name = $self->somatic_variation_build->subject_name;
            $self->status_message("Somatic variation build given. Setting sample name to $sample_name");
            $self->sample_name($sample_name);
        }
    }
    else {
        unless (defined($self->sample_name) && defined($self->input_tsv_file)) {
            die $self->error_message("Sample name and input tsv file must both be defined if no somatic variation build is given")
        }
    }

    unless (-s $self->input_tsv_file) {
        die $self->error_message("Input tsv file %s does not exist or has no size", $self->input_tsv_file);
    }

    unless (Genome::Sys->create_directory($self->output_directory)) {
        die $self->error_message("Coult not create directory (%s)", $self->output_directory);
    }

    # TODO make sure anno db makes sense
    # TODO make sure anno db version makes sense
    # TODO make sure length makes sense
    # TODO make sure allele makes sense
    # TODO make sure epitope_length makes sense
    # TODO make sure netmhc_version makes sense

    return 1;
}

sub _get_workflow_inputs {
    my $self = shift;

    my %inputs = (
        input_tsv_file => $self->input_tsv_file,
        output_directory => $self->output_directory,
        anno_db => $self->anno_db,
        anno_db_version => $self->anno_db_version,
        peptide_sequence_length => $self->peptide_sequence_length,
        allele => $self->allele,
        epitope_length => $self->epitope_length,
        netmhc_version => $self->netmhc_version,
        output_filter => $self->output_filter,
        sample_name => $self->sample_name,
    );

    return \%inputs;
}

sub final_output_file {
    my $self = shift;

    my $file_name = join ('.', $self->sample_name, $self->allele, $self->epitope_length, 'netmhc', 'parsed', $self->output_filter);
    return File::Spec->join($self->output_directory, $file_name);
}

1;
