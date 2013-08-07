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
    has_transient_optional => [
        original_format => { is => 'Text', }, # FIXME should we accept this on the command line?
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Import instrument data...');

    my $original_format = $self->_resolve_original_format;
    return if not $original_format;

    my $method = '_build_workflow_to_import_'.$original_format;
    my $wf = $self->$method;
    return if not $wf;

    my $inputs = $self->_gather_inputs($wf);
    return if not $inputs;

    my $success = Workflow::Simple::run_workflow($wf, %$inputs);
    die 'Run wf failed!' if not $success;

    $self->status_message('Import instrument data...done');
    return 1;
}

sub _resolve_original_format {
    my $self = shift;
    $self->status_message('Resolve original format...');
    
    my @source_files = $self->source_files;
    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my %formats;
    for my $source_file ( @source_files ) {
        my $format = $helpers->source_file_format($source_file);
        return if not $format;
        $formats{$format}++;
    }

    my @formats = keys %formats;
    if ( @formats > 1 ) {
        $self->error_message('Got more than one format when trying to determine format!');
        return;
    }
    $self->status_message('Original format: '.$formats[0]);

    my $max_source_files = ( $formats[0] =~ /^fast[aq]$/ ? 2 : 1 );
    if ( @source_files > $max_source_files ) {
        $self->error_message("Cannot handle more than $max_source_files source files!");
        return;
    }

    $self->status_message('Resolve original format...done');
    return $self->original_format($formats[0]);
}

sub _gather_inputs {
    my ($self, $workflow) = @_;

    Carp::confess('No work flow to gather inputs!') if not $workflow;
    
    my @instrument_data_properties = $self->instrument_data_properties;
    for my $property (qw/ import_source_name description /) {
        my $value = $self->$property;
        next if not defined $value;
        push @instrument_data_properties, $property.'='.$value;
    }

    my @source_files = $self->source_files;
    push @instrument_data_properties, 'original_data_path='.join(',', $self->source_files);
    my $source_path_alias = 'source_'.$self->original_format.'_path';
    my $source_paths = $source_files[0];
    if ( @source_files > 1 ) {
        $source_path_alias .= 's';
        $source_paths = \@source_files;
    }

    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    if ( not $tmp_dir ) {
        $self->error_message('Failed to create tmp dir!');
        return;
    }

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $space_available = $helpers->verify_adequate_disk_space_is_available_for_source_files(
        tmp_dir => $tmp_dir,
        source_files => \@source_files,
    );
    return if not $space_available;

    my %possible_inputs = (
        working_directory => $tmp_dir,
        sample => $self->sample,
        sample_name => $self->sample->name,
        $source_path_alias => $source_paths,
        instrument_data_properties => \@instrument_data_properties,
    );
    return { map { $_ => $possible_inputs{$_} } @{$workflow->operation_type->input_properties} };
}

sub _build_workflow_to_import_fastq {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Inst Data',
        input_properties => [qw/ working_directory source_fastq_paths sample sample_name instrument_data_properties /],
        output_properties => [qw/ instrument_data /],
    );

    my $helper = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my $get_fastqs_op = $helper->add_operation_to_workflow($workflow, 'get fastqs');
    for my $property (qw/ working_directory source_fastq_paths /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $get_fastqs_op,
            right_property => $property,
        );
    }

    my $fastqs_to_bam_op = $helper->add_operation_to_workflow($workflow, 'fastqs to bam');
    for my $property (qw/ working_directory sample_name /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $fastqs_to_bam_op,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $get_fastqs_op,
        left_property => 'fastq_paths',
        right_operation => $fastqs_to_bam_op,
        right_property => 'fastq_paths',
    );

    my $sort_bam_op = $helper->add_operation_to_workflow($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $fastqs_to_bam_op,
        left_property => 'bam_path',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_paths',
    );

    my $create_instdata_and_copy_bam = $helper->add_operation_to_workflow($workflow, 'create instrument data and copy bam');
    for my $property (qw/ sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_paths',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );

    $workflow->add_link(
        left_operation => $create_instdata_and_copy_bam,
        left_property => 'instrument_data',
        right_operation => $workflow->get_output_connector,
        right_property => 'instrument_data',
    );

    return $workflow;
}

sub _build_workflow_to_import_bam {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Inst Data',
        input_properties => [qw/ working_directory source_bam_path sample instrument_data_properties /],
        output_properties => [qw/ instrument_data /],
    );

    my $helper = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my $split_bam_op = $helper->add_operation_to_workflow($workflow, 'split bam by read group');
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'source_bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    my $sort_bam_op = $helper->add_operation_to_workflow($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $split_bam_op,
        left_property => 'read_group_bam_paths',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_paths',
    );
    $sort_bam_op->parallel_by('unsorted_bam_path');

    my $create_instdata_and_copy_bam = $helper->add_operation_to_workflow($workflow, 'create instrument data and copy bam');
    for my $property (qw/ sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_paths',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );
    $create_instdata_and_copy_bam->parallel_by('bam_path');

    $workflow->add_link(
        left_operation => $create_instdata_and_copy_bam,
        left_property => 'instrument_data',
        right_operation => $workflow->get_output_connector,
        right_property => 'instrument_data',
    );

    return $workflow;
}

sub _build_workflow_to_import_sra {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Inst Data',
        input_properties => [qw/ working_directory source_sra_path sample instrument_data_properties /],
        output_properties => [qw/ instrument_data /],
    );

    my $helper = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my $get_sra_op = $helper->add_operation_to_workflow($workflow, 'get sra');
    for my $property (qw/ working_directory source_sra_path /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $get_sra_op,
            right_property => $property,
        );
    }

    my $sra_to_bam_op = $helper->add_operation_to_workflow($workflow, 'sra to bam');
    for my $property (qw/ working_directory sra_path /) {
        $workflow->add_link(
            left_operation => $get_sra_op,
            left_property => $property,
            right_operation => $sra_to_bam_op,
            right_property => $property,
        );
    }

    my $split_bam_op = $helper->add_operation_to_workflow($workflow, 'split bam by read group');
    $workflow->add_link(
        left_operation => $sra_to_bam_op,
        left_property => 'bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    my $sort_bam_op = $helper->add_operation_to_workflow($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $split_bam_op,
        left_property => 'read_group_bam_paths',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_paths',
    );

    my $create_instdata_and_copy_bam = $helper->add_operation_to_workflow($workflow, 'create instrument data and copy bam');
    for my $property (qw/ sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_paths',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );

    $workflow->add_link(
        left_operation => $create_instdata_and_copy_bam,
        left_property => 'instrument_data',
        right_operation => $workflow->get_output_connector,
        right_property => 'instrument_data',
    );

    return $workflow;
}

1;

