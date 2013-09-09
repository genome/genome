package Genome::Model::DeNovoAssembly::Command::PrepareInstrumentData;

use strict;
use warnings;

use Genome;

require File::Temp;

class Genome::Model::DeNovoAssembly::Command::PrepareInstrumentData {
    is => 'Command::V2',
    has_input => [
        build => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1 },
    ],
    has_optional => [
        _input_count => { is => 'Number', default_value => 0, },
        _input_bases => { is => 'Number', default_value => 0, },
        _output_count => { is => 'Number', default_value => 0, },
        _output_bases => { is => 'Number', default_value => 0, },
        _original_base_limit => { is => 'Number', },
        _base_limit => { is => 'Number', },
    ],
    has_constant => [
        lsf_resource => { is => 'Text',
            default_value => "-R 'select[type==LINUX64 && mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]' -M 32000000",
        },
    ],
};

sub _tempdir {
    my $self = shift;

    unless ( $self->{_tempdir} ) {
        $self->{_tempdir} = Genome::Sys->base_temp_directory;
    }

    return $self->{_tempdir};
}

sub _input_metrics_file {
    return $_[0]->_tempdir.'/metrics.input.txt';
}

sub _output_metrics_file {
    return $_[0]->_tempdir.'/metrics.output.txt';
}

sub execute {
    my $self = shift;

    $self->status_message('Prepare instrument data for '.$self->build->description);

    $self->status_message('Verify instrument data...');
    my @instrument_data = $self->build->instrument_data;
    unless ( @instrument_data ) {
        $self->error_message("Failed to prepare instrument data. Build does not have any.");
        return;
    }
    $self->status_message('Verify instrument data...OK');

    $self->_setup_base_limit;

    my @existing_assembler_input_files = $self->build->existing_assembler_input_files;
    if ( @existing_assembler_input_files ) {
        $self->status_message('Removing existing assembler input files');
        for my $file ( @existing_assembler_input_files ) {
            unlink $file;
            if ( -e $file ) {
                $self->error_message("Cannot remove existing assembler input file $file");
                return;
            }
        }
    }

    $self->status_message('Process instrument data');
    INST_DATA: for my $instrument_data ( reverse @instrument_data ) {
        my $process_ok = $self->_process_instrument_data($instrument_data);
        return if not $process_ok;
        my $update_metrics = $self->_update_metrics;
        return if not $update_metrics;
        last INST_DATA if $self->_has_base_limit_been_reached;
    }
    $self->status_message('Process instrument data...OK');

    $self->status_message('Verify assembler input files');
    @existing_assembler_input_files = $self->build->existing_assembler_input_files;
    if ( not @existing_assembler_input_files ) {
        $self->error_message('No assembler input files were created!');
        return;
    }
    $self->status_message('Verify assembler input files...OK');

    my $reads_attempted = $self->_input_count;
    my $reads_processed = $self->_output_count;
    my $reads_processed_success = ( $reads_attempted ? sprintf('%0.3f', $reads_processed / $reads_attempted) : 0);
    $self->build->add_metric(name => 'reads_attempted', value => $reads_attempted);
    $self->build->add_metric(name => 'reads_processed', value => $reads_processed);
    $self->build->add_metric(name => 'reads_processed_success', value => $reads_processed_success);
    $self->status_message('Reads attempted: '.$reads_attempted);
    $self->status_message('Reads processed: '.$reads_processed);
    $self->status_message('Reads processed success: '.($reads_processed_success * 100).'%');

    $self->status_message('Prepare instrument data...OK');
    return 1;
}

sub _setup_base_limit {
    my $self = shift;

    my $base_limit = $self->build->calculate_base_limit_from_coverage;
    return 1 if not defined $base_limit;

    $self->status_message('Setting base limit to: '.$base_limit);
    $self->_original_base_limit($base_limit);
    $self->_base_limit($base_limit);

    return 1;
}

sub _update_metrics {
    my $self = shift;
    $self->status_message('Update metrics...');

    for my $type (qw/ input output /) {
        my $metrics_file_method = '_'.$type.'_metrics_file';
        my $metrics_file = $self->$metrics_file_method;
        $self->status_message(ucfirst($type)." file: $metrics_file");
        if ( not -s $metrics_file ) {
            Carp::confess("No metrics file ($metrics_file) from read processor command.");
        }

        my $metrics = Genome::Model::Tools::Sx::Metrics->from_file($metrics_file);
        if ( not $metrics ) {
            $self->error_message('Failed to create sx metrics from file: '.$metrics_file);
            return;
        }

        for my $name ( $metrics->metric_names ) {
            my $metric_method = '_'.$type.'_'.$name;
            my $metric = $self->$metric_method;
            my $new_metric = $metric + $metrics->$name;
            $self->$metric_method($new_metric);
            $self->status_message("Update $type $name from $metric to $new_metric");
        }
    }

    $self->status_message('Update metrics...OK');
    return 1;
}

sub _has_base_limit_been_reached {
    my $self = shift;

    return if not defined $self->_base_limit;

    $self->status_message('Original base limit: '.$self->_original_base_limit);
    $self->status_message('Bases processed: '.$self->_output_bases);
    my $current_base_limit = $self->_original_base_limit - $self->_output_bases;
    $self->_base_limit($current_base_limit);
    if ( $current_base_limit <= 0 ) {
        $self->status_message('Reached base limit. Stop processing!');
        return 1;
    }
    $self->status_message('New base limit: '.$self->_base_limit);

    $self->status_message('Base limit not reached. Continue processing.');
    return;
}

sub _process_instrument_data {
    my ($self, $instrument_data) = @_;
    $self->status_message('Process: '.join(' ', map { $instrument_data->$_ } (qw/ class id/)));

    # Output files
    my @output_files = $self->build->read_processor_output_files_for_instrument_data($instrument_data);
    return if not @output_files;
    my $output;
    if ( @output_files == 1 ) {
        $output = $output_files[0].':type=sanger:mode=a';
    }
    elsif ( @output_files == 2 ) {
        $output = $output_files[0].':name=fwd:type=sanger:mode=a,'.$output_files[1].':name=rev:type=sanger:mode=a';
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
        my $qual_type = 'sanger'; # imported will be sanger; check solexa
        if ( $instrument_data->can('resolve_quality_converter') ) {
            my $converter = eval{ $instrument_data->resolve_quality_converter };
            $qual_type = 'illumina';
            if ( not $converter ) {
                $self->error_message('No quality converter for instrument data '.$instrument_data->id);
                return;
            }
            elsif ( $converter eq 'sol2sanger' ) {
                $self->error_message('Cannot process old illumina data! Instrument data '.$instrument_data->id);
                return;
            }
            elsif ( $converter eq 'none' ) { # sanger fastqs
                $qual_type = 'sanger';
            }
        }
        my $instrument_data_tempdir = File::Temp::tempdir(CLEANUP => 1);
        if ( not -d $instrument_data_tempdir ) {
            $self->error_message('Failed to make temp directory for instrument data!');
            return;
        }
        my $cmd = "tar -xzf $archive --directory=$instrument_data_tempdir";
        my $tar = Genome::Sys->shellcmd(cmd => $cmd);
        if ( not $tar ) {
            $self->error_message('Failed extract archive for instrument data '.$instrument_data->id);
            return;
        }
        my @input_files = grep { not -d } glob("$instrument_data_tempdir/*");
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
    my $read_processor = $self->build->processing_profile->read_processor;
    my @read_processor_parts = split(/\s+\|\s+/, $read_processor);

    if ( defined $self->_base_limit ) { # coverage limit by bases
        my $current_base_limit = $self->_base_limit;

        # sort reads by qual and use the best reads
        my $read_count = eval { $instrument_data->read_count; };
        my $read_length = eval { $instrument_data->read_length; };
        if ( $read_count and $read_length ) {
            my $inst_data_input_bases = $read_count * $read_length;
            if ( $inst_data_input_bases > $current_base_limit ) {
                $self->status_message('Sorting input reads by average base quality');
                push @read_processor_parts, 'sort by-avg-qual';
            }
        }

        $self->status_message("Limiting bases by base count of $current_base_limit");
        push @read_processor_parts, 'limit by-bases --bases '.$current_base_limit;
    }

    if ( not @read_processor_parts ) { # essentially a copy, but w/ metrics
        @read_processor_parts = ('');
    }

    my @sx_cmd_parts = map { 'gmt sx '.$_ } @read_processor_parts;
    $sx_cmd_parts[0] .= ' --input '.join(',', @inputs);
    $sx_cmd_parts[0] .= ' --input-metrics '.$self->_input_metrics_file;
    $sx_cmd_parts[$#read_processor_parts] .= ' --output '.$output;
    $sx_cmd_parts[$#read_processor_parts] .= ' --output-metrics '.$self->_output_metrics_file;

    # Run
    my $sx_cmd = join(' | ', @sx_cmd_parts);

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $sx_cmd, set_pipefail => 0); };
    if ( not $rv ) {
        $self->error_message('Failed to execute gmt sx command: '.$@);
        return;
    }

    return 1;
}

1;
