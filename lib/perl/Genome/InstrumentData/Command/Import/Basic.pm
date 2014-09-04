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
        library => {
            is => 'Genome::Library',
            doc => 'Library to use.  It must exist.',
        },
        source_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Source files to import. If importing fastqs, put the file containing the forward [read 1] reads first.',
        },
    ],
    has_optional_input => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        description  => {
            is => 'Text',
            doc => 'Description of the data.',
        },
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
        original_format => {
            is => 'Text',
            valid_values => [qw/ bam fastq sra /],
            doc => 'The original format of the source files. Use if the format cannot be determined from the file extension.',
        },
    ],
    has_optional_constant_calculated => {
        helpers => {
            calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ),
        },
    },
    has_optional_transient => [
        _workflow => {},
        _verify_md5_op => {},
        _working_directory => { is => 'Text', },
        _instrument_data_properties => { is => 'Hash', },
        _new_instrument_data => { is => 'Genome::InstrumentData', is_many => 1 },
    ],
};

sub help_detail {
    return <<HELP;
Import sequence files into GMS. All files will be converted to SAM format and stored as BAM.

Source Files 
 Types       Notes
  FASTQ       Can be remote, tar'd and/or gzipped.
  BAM, SAM    Will be split by read group.
  SRA         Aligned and unaligned reads will be dumped. SRAs are known to produce unreliable BAM files.

Instrument Data Properties
 Name and value pairs to add to the instrument data. Separate name and value with an equals (=)
  and name/value pairs with a comma (,).
  
  Example...set flow_cell_id to 'AXXAX' and the index sequence to 'AATTGGCC' on the created instrument
   data entities:

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

    my $workflow = $self->_build_workflow;
    return if not $workflow;

    my $inputs = $self->_gather_inputs_for_workflow;
    return if not $inputs;

    my $success = Workflow::Simple::run_workflow($workflow, %$inputs);
    die 'Run wf failed!' if not $success;

    $self->_new_instrument_data($success->{instrument_data});
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

sub _build_workflow {
    my $self = shift;

    my $workflow = Workflow::Model->create(
        name => 'Import Instrument Data',
        input_properties => [qw/ analysis_project instrument_data_properties library source_paths working_directory /],
        output_properties => [qw/ instrument_data /],
    );
    $self->_workflow($workflow);

    my $retrieve_source_path_op = $self->_add_retrieve_source_path_op_to_workflow;
    return if not $retrieve_source_path_op;

    my $verify_md5_op = $self->_add_verify_md5_op_to_workflow($retrieve_source_path_op);
    return if not $verify_md5_op;
    $self->_verify_md5_op($verify_md5_op);

    my $method = '_build_workflow_to_import_'.$self->original_format;
    my $previous_op = $self->$method;
    return if not $previous_op;

    my $create_instdata_and_copy_bam_op = $self->_add_create_instdata_and_copy_bam_op_to_workflow($previous_op);
    return if not $create_instdata_and_copy_bam_op;

    return $workflow;
}

sub _add_retrieve_source_path_op_to_workflow {
    my $self = shift;

    my $workflow = $self->_workflow;
    my $retrieve_source_path_op = $self->helpers->add_operation_to_workflow_by_name($workflow, 'retrieve source path');
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

    return $retrieve_source_path_op;
}

sub _add_verify_md5_op_to_workflow {
    my ($self, $retrieve_source_path_op) = @_;

    die 'No retrieve source files operation given!' if not $retrieve_source_path_op;

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->helpers->add_operation_to_workflow_by_name($workflow, 'verify md5');
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

    return $verify_md5_op;
}

sub _add_create_instdata_and_copy_bam_op_to_workflow {
    my ($self, $previous_op) = @_;

    die 'No previous op given to _add_create_instdata_and_copy_bam_op_to_workflow!' if not $previous_op;

    my $workflow = $self->_workflow;
    my $create_instdata_and_copy_bam_op = $self->helpers->add_operation_to_workflow_by_name($workflow, 'create instrument data and copy bam');
    for my $property (qw/ analysis_project library instrument_data_properties /) {
        $workflow->add_link(
            left_operation => $workflow->get_input_connector,
            left_property => $property,
            right_operation => $create_instdata_and_copy_bam_op,
            right_property => $property,
        );
    }

    $workflow->add_link(
        left_operation => $previous_op,
        left_property => ( $previous_op->name eq 'sort bam' ) # not ideal...
        ? 'output_bam_path'
        : 'output_bam_paths',
        right_operation => $create_instdata_and_copy_bam_op,
        right_property => 'bam_paths',
    );

    $workflow->add_link(
        left_operation => $self->_verify_md5_op,
        left_property => 'source_md5',
        right_operation => $create_instdata_and_copy_bam_op,
        right_property => 'source_md5s',
    );
    $create_instdata_and_copy_bam_op->parallel_by('bam_path');

    $workflow->add_link(
        left_operation => $create_instdata_and_copy_bam_op,
        left_property => 'instrument_data',
        right_operation => $workflow->get_output_connector,
        right_property => 'instrument_data',
    );

    return $create_instdata_and_copy_bam_op;
}

sub _gather_inputs_for_workflow {
    my $self = shift;

    return {
        analysis_project => $self->analysis_project,
        instrument_data_properties => $self->_instrument_data_properties,
        library => $self->library,
        source_paths => [ $self->source_files ],
        working_directory => $self->_working_directory,
    };
}

sub _build_workflow_to_import_fastq {
    my $self = shift;

    my $helpers = $self->helpers;
    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my %left_op_and_fastqs_property = (
        left_operation => $verify_md5_op,
        left_property => 'source_path',
    );
    my @source_files = $self->source_files;

    # .tar.gz + Archive::Extract->types
    my $archive_command_name = 'archive to fastqs';
    my $archive_command_class_name = $helpers->work_flow_operation_class_from_name($archive_command_name);
    my @archive_types = $archive_command_class_name->types;
    my $is_archived = (
        @source_files == 1
        && grep { $source_files[0] =~ /\Q.$_\E$/ } @archive_types
    );
    if ( $is_archived ) {
        my $archive_to_fastqs_op = $helpers->add_operation_to_workflow_by_name($workflow, $archive_command_name);
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

    my $fastqs_to_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'fastqs to bam');
    for my $property (qw/ working_directory library /) {
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

    my $sort_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $fastqs_to_bam_op,
        left_property => 'output_bam_path',
        right_operation => $sort_bam_op,
        right_property => 'bam_path',
    );

    return $sort_bam_op;
}

sub _build_workflow_to_import_bam {
    my $self = shift;

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my $helpers = $self->helpers;
    my $sanitize_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sanitize bam');
    $workflow->add_link(
        left_operation => $verify_md5_op,
        left_property => 'source_path',
        right_operation => $sanitize_bam_op,
        right_property => 'bam_path',
    );

    my $sort_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $sanitize_bam_op,
        left_property => 'output_bam_path',
        right_operation => $sort_bam_op,
        right_property => 'bam_path',
    );

    my $split_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'split bam by read group');
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'output_bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    return $split_bam_op;
}

sub _build_workflow_to_import_sra {
    my $self = shift;

    my $workflow = $self->_workflow;
    my $verify_md5_op = $self->_verify_md5_op;

    my $helpers = $self->helpers;
    my $sra_to_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sra to bam');
    for my $property_mapping ( [qw/ working_directory working_directory /], [qw/ source_path sra_path /] ) {
        my ($left_property, $right_property) = @$property_mapping;
        $workflow->add_link(
            left_operation => $verify_md5_op,
            left_property => $left_property,
            right_operation => $sra_to_bam_op,
            right_property => $right_property,
        );
    }
    $workflow->add_link(
        left_operation => $workflow->get_input_connector,
        left_property => 'library',
        right_operation => $sra_to_bam_op,
        right_property => 'library',
    );

    my $sanitize_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sanitize bam');
    $workflow->add_link(
        left_operation => $sra_to_bam_op,
        left_property => 'output_bam_path',
        right_operation => $sanitize_bam_op,
        right_property => 'bam_path',
    );

    my $sort_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'sort bam');
    $workflow->add_link(
        left_operation => $sanitize_bam_op,
        left_property => 'output_bam_path',
        right_operation => $sort_bam_op,
        right_property => 'bam_path',
    );

    my $split_bam_op = $helpers->add_operation_to_workflow_by_name($workflow, 'split bam by read group');
    $workflow->add_link(
        left_operation => $sort_bam_op,
        left_property => 'output_bam_path',
        right_operation => $split_bam_op,
        right_property => 'bam_path',
    );

    return $split_bam_op;
}

1;

