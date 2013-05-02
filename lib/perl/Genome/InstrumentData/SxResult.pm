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
            is_many => 1,
        },
    ],
    has_param => [
        read_processor => {
            is => 'Text',
            doc => 'The string describing the read processor operations',
        },
        coverage => {
            is => 'Number',
            doc => 'Desired amount of sequence in output, expressed as a multiple of the estimated genome size',
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
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            calculate => q{
                return Genome::InstrumentData->get([$self->instrument_data_id]);
            }
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

    $self->_prepare_staging_directory;

    my @instrument_data = $self->instrument_data;
    my $num_inst_data = scalar @instrument_data;
    if ($num_inst_data == 1) {
        my $id = $instrument_data[0];
        $self->status_message('Process instrument data '.$id->__display_name__ );
        my $process_ok = $self->_process_instrument_data($id);
        if(not $process_ok) {
            $self->delete;
            return;
        }
        $self->status_message('Process instrument data...OK');
    }
    else {
        my @input_files;
        my $output_file_count = $self->get_output_file_count;
        for (my $input_number=1; $input_number <= $output_file_count; $input_number++) {
            my $new_input_file = join("/", $self->temp_staging_directory, $input_number);
            push @input_files, $new_input_file;
        }
        for my $id (@instrument_data) {
            my %params = (
                instrument_data_id => [$id->id],
                read_processor => $self->read_processor,
                test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
                coverage => undef,
            );
            my @output_file_config = $self->output_file_config;
            if (scalar @output_file_config > 0) {
                $params{output_file_config} = \@output_file_config;
            }
            else {
                $params{output_file_count} = $self->output_file_count;
                $params{output_file_type} = $self->output_file_type;
            }
            my $result = Genome::InstrumentData::SxResult->get_with_lock(
                %params
            );
            unless ($result) {
                $self->error_message("Could not load SX Result for instrument data.  
                    Please create it before using it in a merged result".$id->id);
                $self->delete;
                return;
            }
            #cat the out file(s) onto the master input file(s)
            my @output_files = $result->read_processor_output_file_paths;
            for (my $input_number=0; $input_number < $output_file_count; $input_number++) {
                my $cmd = "cat ".$output_files[$input_number]." >> ".$input_files[$input_number];
                `$cmd`;
            }
        }
        #run the input files through sx for coverage
        #get the estimated genome size
        #be sure it is the same for all instrument data
        my $read_processor = '';
        if ($self->coverage) {
            my $estimated_genome_size;
            for my $id (@instrument_data) {
                my $egs = $id->taxon->estimated_genome_size;
                if ((not defined $estimated_genome_size) or (defined $estimated_genome_size and $egs == $estimated_genome_size)) {
                    $estimated_genome_size = $egs;
                }
                elsif (defined $estimated_genome_size and $egs != $estimated_genome_size) {
                    $self->error_message("Estimated genome size must be the same for all instrument data");
                    $self->delete;
                    return;
                }
            }

            unless (defined $estimated_genome_size) {
                $self->error_message("No estimated genome size found.  This must be specified on the instrument data if a coverage limit is needed for this result");
                $self->delete;
                return;
            }

            #convert $self->coverage and estimated_genome_size to number of bp
            my $bp_count = $self->coverage * $estimated_genome_size;
            #figure out total incoming bases
            #my $incoming_bases = ?;
            #$read_processor = "limit by-bases --select-random-sequences $incoming_bases --bases $bp_count";
            $read_processor = "limit by-bases --bases $bp_count";
        }

        my $output = $self->_output_config;
        unless (defined $output) {
            $self->error_message("Failed to define output");
            $self->delete;
            return;
        }
        my $input_file_type = $self->resolve_merged_input_type;
        my @inputs = map { $_.':type='.$input_file_type } @input_files;
        unless ($self->_run_sx([$read_processor], \@inputs, $output)) {
            $self->error_message("Failed to run sx");
            $self->delete;
            return;
        }
    }

    $self->status_message('Verify output files...');
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
        $self->delete;
        return;
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
    $self->status_message('Verify output files...OK');

    $self->status_message('Process instrument data...OK');

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub resolve_merged_input_type {
    my $self = shift;
    if (defined $self->output_file_type) {
        return $self->output_file_type;
    }
    my @output_configs = $self->output_file_config;
    if (@output_configs) {
       my @parts = split /:/, $output_configs[0];
       for my $part (@parts) {
           if ($part =~ /type=/) {
               $part =~ s/type=//;
               return $part;
           }
       }
    }
    return;
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
    my $base = $self->resolve_base_name_from_instrument_data($self->instrument_data_id);
    return "$base.output_metrics";
}

sub read_processor_input_metric_file {
    my $self = shift;
    my $base = $self->resolve_base_name_from_instrument_data($self->instrument_data_id);
    return "$base.input_metrics";
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

sub _process_instrument_data {
    my ($self, $instrument_data) = @_;
    $self->status_message('Process: '.join(' ', map { $instrument_data->$_ } (qw/ class id/)));

    # Output
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
        my $qual_type = eval{ $instrument_data->native_qual_format; };
        $qual_type //= 'sanger';

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

    return $self->_run_sx(\@read_processor_parts, \@inputs, $output);
}

sub _run_sx {
    my ($self, $read_processor_parts, $inputs, $output) = @_;
    my @sx_cmd_parts = map { 'gmt sx '.$_ } @{$read_processor_parts};
    $sx_cmd_parts[0] .= ' --input '.join(',', @{$inputs});
    $sx_cmd_parts[0] .= ' --input-metrics '.
    $self->temp_staging_directory.'/'.
    $self->read_processor_input_metric_file;
    my $num_parts = scalar @{$read_processor_parts};
    $sx_cmd_parts[$num_parts-1] .= ' --output '.$output;
    $sx_cmd_parts[$num_parts-1] .= ' --output-metrics '.
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

sub resolve_base_name_from_instrument_data {
    my $self = shift;
    my @ids = @_;
    my $base;
    my $num_ids = scalar @ids;
    if ($num_ids == 1) {
        $base = $ids[0];
    }
    else {
        $base = "merged";
        for my $id (@ids) {
            $base .= "-$id";
        }
    }
    return $base;
}

1;

