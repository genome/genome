package Genome::InstrumentData::Command::Import::CompleteGenomics;

use strict;
use warnings;

use Genome;
use Workflow;
use Workflow::Simple;

class Genome::InstrumentData::Command::Import::CompleteGenomics {
    is => 'Command::V2',
    has_input => [
        input_directory => {
            is => 'Text',
            doc => 'path to the Complete Genomics data directory',
        },
        reference_build => { #TODO can this be determined from directories somehow?
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'The reference used to generate the alignment data',
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'The sample that was processed',
        },
    ],
    has_optional_input => [
        verify_integrity => {
            is => 'Boolean',
            default => 1,
            doc => 'compare the sha256sums in manifest.all with the directory contents',
        },
        log_directory => {
            is => 'Text',
            doc => 'Where to write the workflow logs (a temporary working directory is used if not specified)',
        },
        library => {
            is => 'Genome:::Library',
            doc => 'The library of the sample that was processed (will default to a generated library)',
        },
        description => {
            is => 'Text',
            doc => 'A general description of the import data',
        },
        target_region_set_name => {
            is => 'Text',
            doc => 'The target region with which to associate this data (if applicable)',
        },
    ],
    has_transient_optional => [
        _working_directory => {
            is => 'Text',
            doc => 'A network-accessible location to which temporary files can be written (derived from allocation path)',
        },
        instrument_data_id => {
            is => 'Number',
            doc => 'The ID of the newly-created instrument data',
            is_output => 1,
        },
        _instrument_data => {
            is => 'Genome::InstrumentData::Imported',
            id_by => 'instrument_data_id',
            doc => 'The newly-created instrument data',
        },
    ],
    doc => 'A command to bring CompleteGenomics data into the genome modeling system',
};

sub help_detail {
    return 'This command operates on an entire CompleteGenomics directory and converts the sequence data to a BAM file ' .
    'which is then imported as instrument data.  This process takes several hours to complete.';
}

sub execute {
    my $self = shift;

    #setup for the import
    $self->_validate_inputs;
    $self->_create_instrument_data_object;
    $self->_prepare_output_directory;

    #verify the integrity of the directories
    if($self->verify_integrity) {
        return unless $self->_check_directory_integrity;
    }

    #run the individual SAM converters #TODO possibly combine these two workflows for more parallelization?
    $self->status_message('Converting individual files to BAM...');
    $self->_convert_map_to_bam;
    $self->_convert_evidence_to_bam;

    #merge the individual SAM files
    $self->_merge_bam_files;

    #move the finished product to the final location
    $self->_promote_data;

    return 1;
}

sub _validate_inputs {
    my $self = shift;

    $self->status_message('Validating inputs...');

    unless(-d $self->input_directory) {
        die $self->error_message('Could not find directory: ' . $self->input_directory);
    }

    unless(-e $self->input_directory . '/manifest.all') {
        die $self->error_message('Could not find manifest.all in provided directory');
    }

    my $ref = $self->reference_build;
    my $ref_file = $ref->full_consensus_path('crr');
    unless(-e $ref_file) {
        die $self->error_message('Could not find crr file for specified reference.');
    }

    unless($self->library) {
        my $sample = $self->sample;

        my $library = Genome::Library->get(
            name => $sample->name . '-extlibs',
            sample_id => $sample->id
        );
        unless ($library) {
            $library = Genome::Library->create(
                name => $sample->name . '-extlibs',
                sample_id => $sample->id
            );
        }

        unless($library) {
            die $self->error_message('Could not get or create a library for the provided sample.');
        }

        $self->library($library);
    }

    $self->status_message('using library: ' . $self->library->name);
    return 1;
}

sub _create_instrument_data_object {
    my $self = shift;

    my %params = (
        import_format => 'bam',
        import_source_name => 'Complete Genomics',
        is_paired_end => 1,
        library_id => $self->library->id,
        original_data_path => $self->input_directory,
        reference_sequence_build_id => $self->reference_build->id,
        sequencing_platform => 'cg',
    );
    
    if($self->target_region_set_name){
        $params{target_region_set_name} = $self->target_region_set_name;
    }

    my $i = Genome::InstrumentData::Imported->create(%params);
    $self->_instrument_data($i);

    return 1;
}

sub _prepare_output_directory {
    my $self = shift;

    my $instrument_data = $self->_instrument_data;

    my $alloc_path = 'alignment_data/imported/' . $instrument_data->id;
    my $kb_usage = $instrument_data->calculate_alignment_estimated_kb_usage;

    my $alloc = Genome::Disk::Allocation->create(
        disk_group_name     => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
        allocation_path     => $alloc_path,
        kilobytes_requested => ($ENV{UR_DBI_NO_COMMIT}? 5 : $kb_usage),
        owner_class_name    => $instrument_data->class,
        owner_id            => $instrument_data->id,
    );

    my $working_dir = File::Temp::tempdir('working-dir-XXXXX', DIR => $alloc->absolute_path, CLEANUP => 1);
    $self->_working_directory($working_dir);

    return 1;
}

sub _check_directory_integrity {
    my $self = shift;

    $self->status_message('Checking files against manifest file...');

    my $directory = $self->input_directory;
    my $manifest_file = $self->input_directory . '/manifest.all';

    #dies on error
    Genome::Sys->shellcmd(
        cmd => 'cd ' . $directory . '; sha256sum --check ' . $manifest_file,
        input_files => [$manifest_file],
    );

    return 1;
}

sub _generate_conversion_workflow {
    my $self = shift;

    my %options = @_;
    my $name = $options{name};
    my $input_properties = $options{input_properties};
    my $output_properties = $options{output_properties};
    my $operation_command = $options{operation_command};
    my $parallel_by = $options{parallel_by} if defined $options{parallel_by};

    my $wf = Workflow::Model->create(
        name => $name,
        input_properties => $input_properties,
        output_properties => $output_properties,
    );
    $wf->log_dir($self->log_directory);

    my $op = Workflow::Operation->create(
        name => 'converter',
        operation_type => Workflow::OperationType::Command->get($operation_command),
    );
    $op->workflow_model($wf);
    $op->parallel_by($parallel_by);

    for my $prop (@$input_properties) {
        $wf->add_link(
            left_operation => $wf->get_input_connector,
            left_property => $prop,
            right_operation => $op,
            right_property => $prop,
        );
    }

    for my $prop (@$output_properties) {
        $wf->add_link(
            left_operation => $op,
            left_property => $prop,
            right_operation => $wf->get_output_connector,
            right_property => $prop,
        );
    }

    return $wf;
}

sub _convert_map_to_bam {
    my $self = shift;

    #find all the files to convert...
    my @map_file_list = glob(join("/", $self->input_directory, '*', 'MAP', '*', 'mapping_*.tsv.bz2'));

    $self->status_message('Found the following mapping files: ' . join(' ', @map_file_list));

    my @input_properties = ('map_file', 'bam_directory', 'reference_file', 'sort_output');
    my @output_properties = ('bam_file');
    my $wf = $self->_generate_conversion_workflow(
        name => 'map to sam conversion',
        parallel_by => 'map_file',
        input_properties => \@input_properties,
        output_properties => \@output_properties,
        operation_command => 'Genome::Model::Tools::CompleteGenomics::MapToSam',
    );

    my $ref = $self->reference_build;
    my $ref_file = $ref->full_consensus_path('crr');

    my $result = eval {
        Workflow::Simple::run_workflow_lsf(
            $wf,
            map_file => \@map_file_list,
            reference_file => $ref_file,
            bam_directory => $self->_working_directory,
            sort_output => 'coordinate'
        );
    };

    my $err = $@;
    if($err or !$result) {
        die $self->error_message('Failed to generate BAM files from MAP files.' . ($err? $err : ''));
    }

    return 1;
}

sub _convert_evidence_to_bam {
    my $self = shift;

    my @evidence_file_list = glob(join("/", $self->input_directory, '*', 'ASM', 'EVIDENCE', 'evidenceDnbs-*.tsv.bz2'));

    $self->status_message('Found the following evidence files: ' . join(' ', @evidence_file_list));

    my @input_properties = ('evidence_file', 'bam_directory', 'reference_file', 'sort_output');
    my @output_properties = ('bam_file');
    my $wf = $self->_generate_conversion_workflow(
        name => 'evidence to sam conversion',
        parallel_by => 'evidence_file',
        input_properties => \@input_properties,
        output_properties => \@output_properties,
        operation_command => 'Genome::Model::Tools::CompleteGenomics::EvidenceToSam',
    );

    my $ref = $self->reference_build;
    my $ref_file = $ref->full_consensus_path('crr');

    my $result = eval {
        Workflow::Simple::run_workflow_lsf(
            $wf,
            evidence_file => \@evidence_file_list,
            reference_file => $ref_file,
            bam_directory => $self->_working_directory,
            sort_output => 'coordinate'
        );
    };

    my $err = $@;
    if($err or !$result) {
        die $self->error_message('Failed to generate BAM files from Evidence files.' . ($err? $err : ''));
    }

    return 1;
}

sub _merge_bam_files {
    my $self = shift;
    my @bam_file_list = glob(join('/', $self->_working_directory, '*.bam'));
    my $result_bam_file = join('/', $self->_working_directory, 'all_sequences.bam');

    my $tmp_dir = join('/', $self->_working_directory, 'picard-temp');
    Genome::Sys->create_directory($tmp_dir);

    my $merge_cmd = Genome::Model::Tools::Picard::MergeSamFiles->create(
        input_files => \@bam_file_list,
        output_file => $result_bam_file,
        sort_order => 'coordinate',
        temp_directory => $tmp_dir,
        assume_sorted => 1,
    );

    unless($merge_cmd->execute) {
        die $self->error_message('Merge failed.');
    }

    return 1;
}

sub _promote_data {
    my $self = shift;

    my $result_bam_file = join('/', $self->_working_directory, 'all_sequences.bam');
    my $id = $self->_instrument_data;
    my $final_destination = join('/', $id->allocations->absolute_path, 'all_sequences.bam');
    
    unless(rename($result_bam_file, $final_destination)) {
        die $self->error_message('Failed to promote file to final location.');
    }

    return 1;
}

1;
