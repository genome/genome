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
            is_optional => 1,
            valid_values => [1,2],
            doc => 'The number of output files to write to',
        },
        output_file_type => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['sanger','illumina','phred','fasta','bed'],
            doc => 'The type of output file to produce',
        },
        output_file_config => {
            is => 'Text',
            is_many => 1,
            is_optional => 1,
            doc => 'The SX output file config. See "gmt sx --h" for help.',
        },
    ],
    has_metric => [
        input_bases => { is => 'Text', doc => 'Number of bases from input.', },
        input_count => { is => 'Text', doc => 'Number of sequences from input.', },
        output_bases => { is => 'Text', doc => 'Number of bases in output.', },
        output_count => { is => 'Text', doc => 'Number of sequences in output.', },
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

sub _error {
    my ($self, $msg) = @_;
    $self->error_message($msg);
    $self->delete;
    return;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $prepare_staging_directory = $self->_prepare_staging_directory;
    return $self->_error('Failed to prepare tagin directory!') if not $prepare_staging_directory;

    my @sx_command_parts = $self->_construct_sx_command_parts;
    return $self->_error('Failed to run SX!') if not @sx_command_parts; 

    my $run_sx = $self->_run_sx(@sx_command_parts);
    return $self->_error('Failed to run SX!') if not $run_sx; 

    my $set_metrics = $self->set_metrics;
    return $self->_error('Failed to set metrics!') if not $set_metrics; 

    my $verify_output_files = $self->_verify_output_files;
    return $self->_error('Failed to verify SX output files!') if not $verify_output_files; 

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _construct_sx_command_parts {
    my $self = shift;

    my $instrument_data = $self->instrument_data;
    $self->status_message('Instrument data: '.$instrument_data->id);

    # Read processor
    my $read_processor = $self->read_processor;
    $read_processor = '' if not $read_processor;

    # Output files
    my $output = $self->_output_config;
    return if not $output;

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
    elsif ( my $archive = eval{ $instrument_data->attributes(attribute_label => 'archive_path')->attribute_value; } ){
        my $qual_type = eval{ $instrument_data->native_qual_type };
        if ( not $qual_type ) {
            $self->error_message('Failed to get quality type for fastq archive path for instrument data! '.$instrument_data->id);
            return;
        }
        if ( $qual_type eq 'solexa' ) {
            $self->error_message('Cannot process old "solexa" quality type for fastq archive path for instrument data! '.$instrument_data->id);
            return;
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

    return ( $read_processor, \@inputs, $output );
}

sub _run_sx {
    my ($self, $read_processor, $inputs, $output) = @_;

    Carp::confess('No read processor sent to run sx!') if not defined $read_processor;
    Carp::confess('No inputs sent to run sx!') if not $inputs;
    Carp::confess('No output sent to run sx!') if not $output;

    my @read_processor_parts = ( $read_processor ? split(/\s+\|\s+/, $read_processor) : ('') );
    my @sx_cmd_parts = map { 'gmt sx '.$_ } @read_processor_parts;
    $sx_cmd_parts[0] .= ' --input '.join(',', @$inputs);
    $sx_cmd_parts[0] .= ' --input-metrics '.
    $self->temp_staging_input_metric_file;
    my $num_parts = scalar @read_processor_parts;
    $sx_cmd_parts[$num_parts-1] .= ' --output '.$output;
    $sx_cmd_parts[$num_parts-1] .= ' --output-metrics '.
    $self->temp_staging_output_metric_file;

    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->debug_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_handle();
    }

    # Execute command
    my $sx_cmd = join(' | ', @sx_cmd_parts);
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $sx_cmd); };
    if ( not $rv ) {
        $self->error_message('Failed to execute gmt sx command: '.$@);
        return;
    }
    return 1;
}

sub set_metrics {
    my $self = shift;

    my @metric_names = map { $_->property_name } grep { $_->is_metric } $self->__meta__->properties;
    my @metrics_defined = grep { defined $self->$_ } @metric_names;
    return 1 if @metrics_defined == @metric_names;

    $self->status_message('Set metrics...');

    my %metrics = $self->load_metrics;
    return if not %metrics;

    for my $metric_name ( keys %metrics ) {
        $self->$metric_name($metrics{$metric_name});
        $self->status_message( sprintf('%s: %s', ucfirst( join(' ', split('_', $metric_name))), $self->$metric_name) );
    }

    $self->status_message('Set metrics...OK');
    return 1;
}

sub load_metrics {
    my $self = shift;

    my %metrics;
    for my $type (qw/ input output /) {
        my $metric_file_method = 'read_processor_'.$type.'_metric_file';
        my $metric_file = $self->$metric_file_method;
        if ( not $metric_file or not -s $metric_file ) {
            $metric_file_method = 'temp_staging_'.$type.'_metric_file';
            $metric_file = $self->$metric_file_method;
            if ( not -s $metric_file ) {
                $self->error_message(ucfirst($type).' metric file not created!');
                return;
            }
        }

        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->from_file($metric_file);
        if ( not $metrics ) {
            $self->error_message('Failed to create SX metrics from file! '.$metric_file);
            return;
        }

        for my $metric_name (qw/ bases count /) {
            my $metric_value = $metrics->$metric_name;
            if ( not defined $metric_value ) {
                $self->error_message("No $metric_name in metrics! ".Data::Dumper::Dumper($metrics));
                return;
            }
            $metrics{$type.'_'.$metric_name} = $metric_value;
        }
    }

    return %metrics;
}

sub _verify_output_files {
    my $self = shift;
    $self->status_message('Verify output files...');

    my $output_count = $self->output_count;
    if ( $output_count == 0 ) {
        $self->status_message('SX ran correctly, but all sequences were filtered out. The output files are empty.');
        return 1;
    }

    # TODO run a  validator?
    my @output_files = $self->read_processor_output_files;
    my $existing_cnt = 0;
    foreach my $output_file (@output_files) {
        my $temp_staging_output_file = $self->temp_staging_directory.'/'.$output_file;
        if ( not -s $temp_staging_output_file ) {
            unlink $temp_staging_output_file;
        }
        else {
            $existing_cnt++;
        }
    }
    if ( not $existing_cnt ) {
        $self->error_message('No output files were created!');
        return;
    }

    $self->status_message('Verify output files...OK');
    return 1;
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
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub resolve_base_name_from_instrument_data {
    my $self = shift;
    return $self->instrument_data_id;
}

sub read_processor_input_metric_file_base_name {
    my $self = shift;
    my $base = $self->resolve_base_name_from_instrument_data;
    return "$base.input_metrics";
}

sub temp_staging_input_metric_file {
    my $self = shift;
    my $temp_staging_directory = $self->temp_staging_directory;
    return if not $temp_staging_directory or not -d $temp_staging_directory;
    return $temp_staging_directory.'/'.$self->read_processor_input_metric_file_base_name;
}

sub read_processor_input_metric_file {
    my $self = shift;
    my $output_dir = $self->output_dir;
    return if not $output_dir or not -d $output_dir;
    return $output_dir.'/'.$self->read_processor_input_metric_file_base_name;
}

sub read_processor_output_metric_file_base_name {
    my $self = shift;
    my $base = $self->resolve_base_name_from_instrument_data;
    return "$base.output_metrics";
}

sub temp_staging_output_metric_file {
    my $self = shift;
    my $temp_staging_directory = $self->temp_staging_directory;
    return if not $temp_staging_directory or not -d $temp_staging_directory;
    return $temp_staging_directory.'/'.$self->read_processor_output_metric_file_base_name;
}

sub read_processor_output_metric_file {
    my $self = shift;
    my $output_dir = $self->output_dir;
    return if not $output_dir or not -d $output_dir;
    return $output_dir.'/'.$self->read_processor_output_metric_file_base_name;
}

sub read_processor_output_files {
    my $self = shift;

    my @output_config_params = $self->_output_config_params;
    return if not @output_config_params;

    return ( map { $_->{basename} } @output_config_params );
}

sub read_processor_output_file_paths {
    my $self = shift;

    my @read_processor_output_files = $self->read_processor_output_files;
    return if not @read_processor_output_files;

    return ( map { $self->output_dir.'/'.$_ } @read_processor_output_files );
}


sub _output_config {
    my $self = shift;

    my @output_config_params = $self->_output_config_params;
    return if not @output_config_params;

    for my $params ( @output_config_params ) {
        $params->{file} = $self->temp_staging_directory.'/'.delete($params->{basename});
        for my $key ( keys %$params ) {
            next if $key !~ /_file$/;
            $params->{$key} = $self->temp_staging_directory.'/'.$params->{$key};
        }
    }

    return join(',', map { Genome::Model::Tools::Sx::Functions->hash_to_config(%$_) } @output_config_params);
}

sub _output_config_params {
    my $self = shift;

    my @params;
    if ( my @output_file_config = $self->output_file_config ) {
        for my $config ( @output_file_config ) {
            my %params = Genome::Model::Tools::Sx::Functions->config_to_hash($config);
            if ( not $params{basename} ) {
                $self->error_message('No basename for SX output config params! '.Data::Dumper::Dumper(\%params));
                return;
            }
            elsif ( $params{basename} !~ /^[\w\d\.\-_]+$/ ) {
                $self->error_message('Malformed basename for SX output config params! '.Data::Dumper::Dumper(\%params));
                return;
            }
            if ( not $params{type} ) {
                $self->error_message('No type for SX output config params! '.Data::Dumper::Dumper(\%params));
                return;
            }
            push @params, \%params;
        }
    }
    elsif ( $self->output_file_count ) {
        my $base = $self->resolve_base_name_from_instrument_data($self->instrument_data_id);

        if ( $self->output_file_count == 1 ) {
            @params = ({
                    basename => $base.'.'.$self->output_file_suffix,
                    type => $self->output_file_type,
                });
        }
        elsif ( $self->output_file_count == 2 ) {
            @params = ({
                    basename => $base.'.1.'.$self->output_file_suffix,
                    name => 'fwd',
                    type => $self->output_file_type,
                },
                {
                    basename => $base.'.2.'.  $self->output_file_suffix,
                    name => 'rev',
                    type => $self->output_file_type,
                });
        }
    }
    else {
        $self->error_message('No output config or output file count set!');
        return;
    }

    for my $params ( @params ) {
        $params->{mode} = 'w';
    }

    return @params;
}

sub get_output_file_count {
    my $self = shift;
    my @configs = $self->output_file_config;
    if ($self->output_file_count) {
        return $self->output_file_count;
    }
    elsif (scalar @configs >= 1) {
        return scalar @configs;
    }
    else {
        $self->error_message("Couldn't resolve output file count");
        die $self->error_message;
    }
}

1;

