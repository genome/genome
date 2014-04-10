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

sub execute {
    my $self = shift;

    $self->debug_message("Validating Inputs...");
    $self->_validate_inputs();

    $self->debug_message("Constructing Workflow...");
    my $workflow = $self->_construct_workflow();

    $self->debug_message("Getting Workflow Inputs...");
    my $inputs = $self->_get_workflow_inputs();

    $self->debug_message("Running Workflow...");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %$inputs);

    unless($result){
        $self->error_message( join("\n", map($_->name . ': ' . $_->error, @Workflow::Simple::ERROR)) );
        die $self->error_message("Workflow did not return correctly.");
    }

    return 1;
}

sub _construct_workflow {
    my ($self) = @_;

    my $xml = __FILE__ . '.xml';
    my $workflow = Workflow::Operation->create_from_xml($xml);
    $workflow->log_dir($self->output_directory);

    return $workflow;
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
        length => $self->peptide_sequence_length,
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
