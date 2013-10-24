package Genome::InstrumentData::Command::Import::Basic;

use strict;
use warnings;

use Genome;

use Workflow::Simple;

class Genome::InstrumentData::Command::Import::Basic { 
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
    has_optional_param => [
        original_format => {
            is => 'Text',
            valid_values => [qw/ bam fastq sra /],
            doc => 'The original format of the source files. Use if the format cannot be determined from the file extension.',
        },
    ],
    has_optional_transient => [
        _workflow => {},
        _verify_md5_op => {},
        _working_directory => { is => 'Text', },
        _instrument_data_properties => { is => 'Hash', },
    ],
};

sub help_detail {
    return <<HELP;
    Import sequence files into GMS. All files will be converted to SAM format and stored as BAM.


    Source Files 

     Types      Notes
     FASTQ      Can be remote, tar'd and/or gzipped.
     SAM [BAM]  Will be split by read group.
     SRA        Aligned and unaligned reads will be dumped. SRAs are known to produce unreliable BAM files.


     Instrument Data Properties
      Indicate properties for the resulting instrument data entity. Give as comma separtated key=values pairs.

      Examples:
       flow_cell_id=AXXAX,index_sequence=AATTGGCC

HELP
}

sub execute {
    my $self = shift;
    $self->status_message('Import instrument data...');

    my $instdat_props_ok = $self->_resolve_instrument_data_properties;
    return if not $instdat_props_ok;

    my $original_format = $self->_resolve_original_format;
    return if not $original_format;

    my $working_directory = $self->_resolve_working_directory;
    return if not $working_directory;

    my $space_available = $self->_verify_adequate_disk_space_is_available_for_source_files;
    return if not $space_available;

    my $workflow = $self->_create_workflow;
    return if not $workflow;

    my $method = '_build_workflow_to_import_'.$original_format;
    my $wf = $self->$method;
    return if not $wf;

    my $inputs = $self->_gather_inputs_for_workflow;
    return if not $inputs;

    my $success = Workflow::Simple::run_workflow($wf, %$inputs);
    die 'Run wf failed!' if not $success;

    $self->status_message('Import instrument data...done');
    return 1;
}

sub _resolve_original_format {
    my $self = shift;
    $self->status_message('Resolve original format...');

    if ( $self->original_format ) {
        $self->status_message('Original format: '.$self->original_format);
        return $self->original_format;
    }

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

sub _resolve_instrument_data_properties {
    my $self = shift;

    my $properties = {};
    if ( $self->instrument_data_properties ) {
        my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
        $properties = $helpers->key_value_pairs_to_hash( $self->instrument_data_properties );
        return if not $properties;
    }

    for my $name (qw/ import_source_name description /) {
        my $value = $self->$name;
        next if not defined $value;
        $properties->{$name} = $value;
    }

    $properties->{original_data_path} = join(',', $self->source_files);

    return $self->_instrument_data_properties($properties);
}

sub _resolve_working_directory {
    my $self = shift;

    my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
    if ( not $tmp_dir ) {
        $self->error_message('Failed to create tmp dir!');
        return;
    }

    return $self->_working_directory($tmp_dir);
}

sub _verify_adequate_disk_space_is_available_for_source_files {
    my $self = shift;
    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $space_available = $helpers->verify_adequate_disk_space_is_available_for_source_files(
        working_directory => $self->_working_directory,
        source_files => [ $self->source_files ],
    );
    return $space_available;
}

sub _create_workflow {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Instrument Data',
        input_properties => [qw/ working_directory source_paths sample instrument_data_properties /],
        output_properties => [qw/ instrument_data /],
    );
    $self->_workflow($workflow);

    my $retrieve_source_path_op = $self->_add_operation_to_workflow('retrieve source path');
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'working_directory',
        right_operation => $retrieve_source_path_op,
        right_property => 'working_directory',
    );
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'source_paths',
        right_operation => $retrieve_source_path_op,
        right_property => 'source_path',
    );
    $retrieve_source_path_op->parallel_by('source_path') if $self->source_files > 1;

    my $verify_md5_op = $self->_add_operation_to_workflow('verify md5');
    $self->_verify_md5_op($verify_md5_op);
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'working_directory',
        right_operation => $verify_md5_op,
        right_property => 'working_directory',
    );
    $workflow->add_link(
        left_operation => $retrieve_source_path_op,
        left_property => 'destination_path',
        right_operation => $verify_md5_op,
        right_property => 'source_path',
   );
   $verify_md5_op->parallel_by('source_path') if $self->source_files > 1;

    return $workflow;
}

sub _add_operation_to_workflow {
    my ($self, $name) = @_;

    my $command_class_name = 'Genome::InstrumentData::Command::Import::WorkFlow::'.join('', map { ucfirst } split(' ', $name));
    my $operation_type = Workflow::OperationType::Command->create(command_class_name => $command_class_name);
    if ( not $operation_type ) {
        $self->error_message("Failed to create work flow operation for $name");
        return;
    }

    my $operation = $self->_workflow->add_operation(
        name => $name,
        operation_type => $operation_type,
    );

    return $operation;
}

sub _gather_inputs_for_workflow {
    my $self = shift;

    return {
        working_directory => $self->_working_directory,
        sample => $self->sample,
        source_paths => [ $self->source_files ],
        instrument_data_properties => $self->_instrument_data_properties,
    };
}

sub _build_workflow_to_import_fastq {
    my $self = shift;

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my %left_op_and_fastqs_property = (
        left_operation => $verify_md5_op,
        left_property => 'source_path',
    );
    my @source_files = $self->source_files;
    if ( @source_files == 1 and $source_files[0] =~ /\.t?gz$/ ){
        my $archive_to_fastqs_op = $self->_add_operation_to_workflow('archive to fastqs');
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => 'working_directory',
            right_operation => $archive_to_fastqs_op,
            right_property => 'working_directory',
        );
        $workflow->add_link(
            left_operation => $verify_md5_op,
            left_property => 'source_path',
            right_operation => $archive_to_fastqs_op,
            right_property => 'archive_path',
        );
        %left_op_and_fastqs_property = (
            left_operation => $archive_to_fastqs_op,
            left_property => 'fastq_paths',
        );
    }

    my $fastqs_to_bam_op = $self->_add_operation_to_workflow('fastqs to bam');
    for my $property (qw/ working_directory sample /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $fastqs_to_bam_op,
            right_property => $property,
        );
    }
    $workflow->add_link(
        %left_op_and_fastqs_property,
        right_operation => $fastqs_to_bam_op,
        right_property => 'fastq_paths',
    );

    my $sort_bam_op = $self->_add_operation_to_workflow('sort bam');
    $workflow->add_link(
        left_operation => $fastqs_to_bam_op,
        left_property => 'bam_path',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_path',
    );

    my $create_instdata_and_copy_bam = $self->_add_operation_to_workflow('create instrument data and copy bam');
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
        left_property => 'sorted_bam_path',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );
    $workflow->add_link(
        left_operation => $verify_md5_op,
        left_property => 'source_md5',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'source_md5s',
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

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my $sort_bam_op = $self->_add_operation_to_workflow('sort bam');
    $workflow->add_link(
        left_operation => $verify_md5_op,
        left_property => 'source_path',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_path',
    );

    my $split_bam_op = $self->_add_operation_to_workflow('split bam by read group');
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    my $create_instdata_and_copy_bam = $self->_add_operation_to_workflow('create instrument data and copy bam');
    for my $property (qw/ sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $split_bam_op,
        left_property => 'read_group_bam_paths',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );
    $workflow->add_link(
        left_operation => $verify_md5_op,
        left_property => 'source_md5',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'source_md5s',
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

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my $sra_to_bam_op = $self->_add_operation_to_workflow('sra to bam');
    for my $property_mapping ( [qw/ working_directory working_directory /], [qw/ source_path sra_path /] ) {
        my ($left_property, $right_property) = @$property_mapping;
        $workflow->add_link(
            left_operation => $verify_md5_op,
            left_property => $left_property,
            right_operation => $sra_to_bam_op,
            right_property => $right_property,
        );
    }

    my $sort_bam_op = $self->_add_operation_to_workflow('sort bam');
    $workflow->add_link(
        left_operation => $sra_to_bam_op,
        left_property => 'bam_path',
        right_operation => $sort_bam_op,
        right_property => 'unsorted_bam_path',
    );

    my $split_bam_op = $self->_add_operation_to_workflow('split bam by read group');
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'sorted_bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    my $create_instdata_and_copy_bam = $self->_add_operation_to_workflow('create instrument data and copy bam');
    for my $property (qw/ sample instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam,
            right_property => $property,
        );
    }
    $workflow->add_link(
        left_operation => $split_bam_op,
        left_property => 'read_group_bam_paths',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'bam_paths',
    );
    $workflow->add_link(
        left_operation => $verify_md5_op,
        left_property => 'source_md5',
        right_operation => $create_instdata_and_copy_bam,
        right_property => 'source_md5s',
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

1;

