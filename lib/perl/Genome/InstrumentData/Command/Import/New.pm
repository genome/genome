package Genome::InstrumentData::Command::Import::New;

use strict;
use warnings;

use Genome;

use Workflow::Simple;

class Genome::InstrumentData::Command::Import::New { 
    is => 'Command::V2',
    has_input => [
        import_source_name => {
            is => 'Text',
            doc => 'Organiztion name or abbreviation from where the source file(s) were generated or downloaded.',
        },
        source_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Source files to import. If importing multiple files, put the file containing the forward reads first.',
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'Sample to use. The external library for the instrument data will be gotten or created.',
        },
    ],
    has_optional_input => [
        description  => {
            is => 'Text',
            doc => 'Description of the data.',
        },
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
    ],
    has_output => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The instrument data entity to be imported.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Import instrument data...');

    my $inputs = $self->_gather_inputs;
    return if not $inputs;

    my $wf = $self->_build_workflow;
    return if not $wf;

    my $success = Workflow::Simple::run_workflow($wf, %$inputs);
    die 'Run wf failed!' if not $success;

    $self->status_message('Import instrument data...done');
    return 1;
}

sub _gather_inputs {
    my $self = shift;
    my @instrument_data_properties = $self->instrument_data_properties;
    for my $property (qw/ import_source_name description /) {
        my $value = $self->$property;
        next if not defined $value;
        push @instrument_data_properties, $property.'='.$value;
    }
    return {
        source_files => [ $self->source_files ],
        sample => $self->sample,
        instrument_data_properties => \@instrument_data_properties,
    };
}

sub _build_workflow {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Inst Data',
        input_properties => [qw/ source_files sample instrument_data_properties /],
        output_properties => [qw/ instrument_data /],
    );

    my $helper = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $create_instrument_data_op = $helper->add_operation_to_workflow($workflow, 'create instrument data');
    for my $property (qw/ source_files sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instrument_data_op,
            right_property => $property,
        );
    }

    my $transfer_fastqs_op = $helper->add_operation_to_workflow($workflow, 'transfer fastqs');
    $workflow->add_link(
        left_operation => $create_instrument_data_op,
        left_property => 'instrument_data',
        right_operation => $transfer_fastqs_op,
        right_property => 'instrument_data',
    );

    my $convert_to_bam_op = $helper->add_operation_to_workflow($workflow, 'convert fastqs to bam');
    $workflow->add_link(
        left_operation => $transfer_fastqs_op,
        left_property => 'instrument_data',
        right_operation => $convert_to_bam_op,
        right_property => 'instrument_data',
    );

    my $sort_bam_op = $helper->add_operation_to_workflow($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $convert_to_bam_op,
        left_property => 'instrument_data',
        right_operation => $sort_bam_op,
        right_property => 'instrument_data',
    );
    $workflow->add_link(
        left_operation => $convert_to_bam_op,
        left_property => 'bam_path',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_path',
    );

    my $verify_and_move_bam_op = $helper->add_operation_to_workflow($workflow, 'verify and move bam to final location');
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'instrument_data',
        right_operation => $verify_and_move_bam_op,
        right_property => 'instrument_data',
    );
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_path',
        right_operation => $verify_and_move_bam_op,
        right_property => 'bam_path',
    );

    $workflow->add_link(
        left_operation => $verify_and_move_bam_op,
        left_property => 'instrument_data',
        right_operation => $workflow->get_output_connector,
        right_property => 'instrument_data',
    );

    return $workflow;
}

1;

