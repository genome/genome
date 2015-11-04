package Genome::InstrumentData::Command::Import::Launch;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::CsvParser;
use Genome::InstrumentData::Command::Import::WorkFlow::Inputs;
use Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles;
require List::Util;

class Genome::InstrumentData::Command::Import::Launch {
    is => 'Command::V2',
    doc => 'Manage importing sequence files into GMS',
    has => {
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        file => {
            is => 'Text',
            doc => 'The metadata file containing sequence files, library names, and other information to be associated with instrument data.',
        },
        job_group_name => {
            is => 'Text',
            doc => 'The job group name. Used to throttle imports to prvent too many running at a time.',
        },
    },
    has_optional => {
        mem => {
            is => 'Number',
            default_value => 8000,
            doc => 'Amount of memory in megabytes to request for each import.',
        },
    },
    has_optional_transient => {
        _wf_inputs => { is => 'Array', },
        gtmp => { is => 'Number', },
        process => { is => 'Genome::InstrumentData::Command::Import::Process', },
    },
};

sub help_brief {
    return 'batch import sequence files into GMS'
}

sub help_detail {
    my $help = <<HELP;
Given a metadata file, launch an import process that will import the sequence files in GMS. The launching of the jobs is handled a genome 'process'.

Listing status of a process:

\$ genome process view \$PROCESS_ID

Listing created instrument data:

\$ genome instrument-data list imported process_id=\$PROCESS_ID

About the Metadata File

HELP
    $help .= Genome::InstrumentData::Command::Import::CsvParser->csv_help;
    return $help;
}

sub execute {
    my $self = shift;

    $self->_check_for_running_processes;
    $self->process( Genome::InstrumentData::Command::Import::Process->create(import_file => $self->file) );
    $self->_create_wf_inputs;
    $self->_calculate_gtmp_required;
    $self->_launch_process;

    return 1
}

sub _check_for_running_processes {
    my $self = shift;

    my $md5 = Genome::Sys->md5sum($self->file);
    die $self->error_message('Failed to get md5 for import file! %s', $self->file) if not $md5;

    my @active_processes = Genome::InstrumentData::Command::Import::Process->get(
        import_md5 => $md5,
        status => [qw/ New Scheduled Running /],
    );

    return 1 if not @active_processes;

    $self->status_message("Found '%s' process (%s) for metadata file: %s", $active_processes[0]->status, $active_processes[0]->id, $self->file);
    die $self->error_message('Cannot start another import process until the previous one has completed!');
}

sub _create_wf_inputs {
    my $self = shift;

    my @inputs;
    my $parser = Genome::InstrumentData::Command::Import::CsvParser->create(file => $self->file);
    my %seen;
    while ( my $import = $parser->next ) {
        my $library_name = $import->{library}->{name};
        my $string = join(' ', $library_name, join(',', $import->{source_files}), map { $import->{instdata}->{$_} } keys %{$import->{instdata}});
        my $id = substr(Genome::Sys->md5sum_data($string), 0, 6);
        if ( $seen{$id} ) {
            die $self->error_message("Duplicate source file/library combination! $string");
        }
        $seen{$id}++;

        my @libraries = Genome::Library->get(name => $library_name);
        die $self->error_message('No library for name: %s', $library_name) if not @libraries;
        die $self->error_message('Multiple libraries for library name: %s', $library_name) if @libraries > 1;

        push @inputs, Genome::InstrumentData::Command::Import::WorkFlow::Inputs->create(
            process_id => $self->process->id,
            analysis_project_id => $self->analysis_project->id,
            library_id => $libraries[0]->id,
            instrument_data_properties => $import->{instdata},
            source_paths => $import->{source_files},
        );
    }

    return $self->_wf_inputs(\@inputs);
}

sub _calculate_gtmp_required {
    my $self = shift;

    my @kb_required;
    for my $input ( @{$self->_wf_inputs} ) {
        my $kb_required = $input->source_files->kilobytes_required_for_processing;
        $kb_required = 1048576 if $kb_required < 1048576; # 1 Gb 
        push @kb_required, $kb_required;
    }

    my $max_kb_required = List::Util::max(@kb_required);
    return $self->gtmp( $max_kb_required / ( 1024 * 1024 ) );
}

sub _launch_process {
    my $self = shift;

    my $dag = Genome::WorkflowBuilder::DAG->create(name => 'Import Instrument Data for '.$self->file);
    my $lsf_resource = sprintf(
        "-g %s -M %s -R 'select [mem>%s & gtmp>%s] rsuage[mem=%s,gtmp=%s]", 
        $self->job_group_name, ($self->mem * 1024), $self->mem, $self->gtmp, $self->mem, $self->gtmp,
    );
    my $import_op = Genome::WorkflowBuilder::Command->create(
        name => 'InstData Import : Run WF',
        command => 'Genome::InstrumentData::Command::Import::WorkFlow::Run',
        lsf_resource => $lsf_resource,
    );
    $dag->connect_input(
        input_property => 'work_flow_inputs',
        destination => $import_op,
        destination_property => 'work_flow_inputs',
    );
    $dag->add_operation($import_op);
    $dag->parallel_by('work_flow_inputs');
    $dag->connect_output(
        output_property => 'instrument_data',
        source => $import_op,
        source_property => 'instrument_data',
    );

    $self->process->run(
        workflow_xml => $dag->get_xml,
        workflow_inputs => { work_flow_inputs => $self->_wf_inputs, },
    );
    $self->status_message("Started imports!\nProcess id: %s\nMetadata directory: %s\nView status with:'genome process view %s'", $self->process->id, $self->process->metadata_directory, $self->process->id);

    return 1;
}

1;

