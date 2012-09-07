package Genome::InstrumentData::SxResult;

use strict;
use warnings;

use Sys::Hostname;
use Genome;

class Genome::InstrumentData::SxResult {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        instrument_data_id => {
            is => 'Text',
            doc => 'The local database id of the instrument data to operate on',
        },
    ],
    has_param => [
        read_processor => {
            is => 'Text',
            doc => 'The string describing the read processor operations',
        },
        output_file_count => {
            is => 'Number',
            valid_values => [1,2],
            doc => 'The number of output files to write to',
        },
        output_file_type => {
            is => 'Text',
            valid_values => ['sanger','illumina','phred','fasta','bed'],
            doc => 'The type of output file to produce',
        },
    ],
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
        },
        output_file_suffix => {
            is => 'Text',
            is_calculated => 1,
            calculate_from => ['output_file_type'],
            calculate => q { 
                if ($output_file_type eq 'sanger' or $output_file_type eq 'illumina' or $output_file_type eq 'phred'){
                    return 'fastq'}
                elsif ($output_file_type eq 'fasta') {
                    return 'fasta'}
                elsif ($output_file_type eq 'bed') {
                    return 'bed'}
            },
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    my $instrument_data = $self->instrument_data;


    $self->_prepare_staging_directory;

    $self->status_message('Process instrument data '.$instrument_data->__display_name__ );

    my $process_ok = $self->_process_instrument_data;
    if(not $process_ok) {
        $self->delete;
        return;
    }

    $self->status_message('Process instrument data...OK');

    $self->status_message('Verify assembler input files');
    my @output_files = $self->read_processor_output_files;
    foreach my $output_file (@output_files) {
        if ( not -s $self->temp_staging_directory.'/'.
                $output_file) {
            $self->error_message('Output file '.$output_file.
                ' was not created');
            $self->delete;
            return;
        }
    }
    if (not -s $self->temp_staging_directory.'/'.
            $self->read_processor_output_metric_file) {
        $self->error_message('Output metrics file not created');
        $self->delete;
        return;
    }
    if (not -s $self->temp_staging_directory.'/'.
            $self->read_processor_input_metric_file) {
        $self->error_message('Input metrics file not created');
        $self->delete;
        return;
    }
    $self->status_message('Verify assembler input files...OK');

    $self->status_message('Process instrument data...OK');

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $hostname = hostname;

    my $user = $ENV{'USER'};
    my $base_dir = sprintf("sxresult-%s-%s-%s-%s",           $hostname, $user, $$, $self->id);
    my $directory = join('/', 'build_merged_alignments',$self->id,$base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    return 'info_genome_models';
}

sub read_processor_output_metric_file {
    my $self = shift;
    return $self->instrument_data_id.'.output_metrics';
}

sub read_processor_input_metric_file {
    my $self = shift;
    return $self->instrument_data_id.'.input_metrics';
}

sub read_processor_output_files {
    my $self = shift;

    my $output_file_count = $self->output_file_count;

    if ($output_file_count == 1) {
        my $file_name = 
            $self->instrument_data_id.'.'.
            $self->output_file_suffix;
        return ($file_name);
    }
    elsif ($output_file_count == 2) {
        my $file_name_1 = 
            $self->instrument_data_id.'.1.'.
            $self->output_file_suffix;
        my $file_name_2 = 
            $self->instrument_data_id.'.2.'.
            $self->output_file_suffix;
        return ($file_name_1,$file_name_2);
    }
    else {
        $self->error_message('Can only handle 1 or 2 output files');
        die $self->error_message();
    }
}

sub _process_instrument_data {
    my ($self) = @_;
    my $instrument_data = $self->instrument_data;
    $self->status_message('Process: '.join(' ', map { $instrument_data->$_ } (qw/ class id/)));

    # Output files
    my $output;
    my $output_file_count = $self->output_file_count;
    my $output_file_type = $self->output_file_type;
    my $output_file_mode = 'w';

    my @output_files = $self->read_processor_output_files;

    if ( $output_file_count == 1 ) {
        $output = $self->temp_staging_directory.'/'.
        $output_files[0].':type='.$output_file_type.':mode='.$output_file_mode;
    }
    elsif ( $output_file_count == 2 ) {
        $output = $self->temp_staging_directory.'/'.
        $output_files[0].':name=fwd:type='.
            $output_file_type.':mode='.$output_file_mode.','
            .$self->temp_staging_directory.'/'.
            $output_files[1].
            ':name=rev:type='.$output_file_type.':mode='.
            $output_file_mode;
    }
    else {
        $self->error_message('Cannot handle more than 2 output files');
        return;
    }

    # Input files
    my $is_paired_end = eval{ $instrument_data->is_paired_end; };
    my $input_cnt = ( $is_paired_end ? 2 : 1 );
    my @inputs;
    if ( my $bam = eval{ $instrument_data->bam_path } ) {
        @inputs = ( $bam.':type=bam:cnt='.$input_cnt );
    }
    elsif ( my $sff = eval{ $instrument_data->sff_file } ) {
        @inputs = ( $sff.':type=sff:cnt='.$input_cnt );
    }
    elsif ( my $archive = eval{ $instrument_data->archive_path; } ){
        my $qual_type = $instrument_data->native_qual_format;
        
        if ( $instrument_data->can('resolve_quality_converter') ) {
            my $converter = eval{ $instrument_data->resolve_quality_converter };
            if ( not $converter ) {
                $self->error_message('No quality converter for instrument data '.$instrument_data->id);
                return;
            }
            elsif ( $converter eq 'sol2sanger' ) {
                $self->error_message('Cannot process old illumina data! Instrument data '.$instrument_data->id);
                return;
            }
        }
        my $cmd = "tar -xzf $archive --directory=".$self->temp_staging_directory;
        my $tar = Genome::Sys->shellcmd(cmd => $cmd);
        if ( not $tar ) {
            $self->error_message('Failed extract archive for instrument data '.$instrument_data->id);
            return;
        }
        my @input_files = grep { not -d } glob($self->temp_staging_directory."/*");
        if ( not @input_files ) {
            $self->error_message('No fastqs from archive from instrument data '.$instrument_data->id);
            return;
        }
        @inputs = map { $_.':type='.$qual_type } @input_files;
    }
    else {
        $self->error_message('Failed to get bam, sff or archived fastqs from instrument data: '.$instrument_data->id);
        return;
    }

    # Sx read processor
    my @read_processor_parts = split(/\s+\|\s+/, $self->read_processor);

    if ( not @read_processor_parts ) { # essentially a copy, but w/ metrics
        @read_processor_parts = ('');
    }

    my @sx_cmd_parts = map { 'gmt sx '.$_ } @read_processor_parts;
    $sx_cmd_parts[0] .= ' --input '.join(',', @inputs);
    $sx_cmd_parts[0] .= ' --input-metrics '.
    $self->temp_staging_directory.'/'.
    $self->read_processor_input_metric_file;
    $sx_cmd_parts[$#read_processor_parts] .= ' --output '.$output;
    $sx_cmd_parts[$#read_processor_parts] .= ' --output-metrics '.
    $self->temp_staging_directory.'/'.
    $self->read_processor_output_metric_file;

    # Run
    my $sx_cmd = join(' | ', @sx_cmd_parts);
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $sx_cmd); };
    if ( not $rv ) {
        $self->error_message('Failed to execute gmt sx command: '.$@);
        return;
    }

    return 1;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id=>$params_bx->id,
        inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

1;

