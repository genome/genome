package Genome::InstrumentData::MergedSxResult;

use strict;
use warnings;
use Genome;

class Genome::InstrumentData::MergedSxResult {
    is => 'Genome::InstrumentData::SxResult',
    has_input => [
        instrument_data_id => {
            is => 'Text',
            doc => 'The local database id of the instrument data to operate on',
            is_many => 1,
        },
    ],
    has_param => [
        coverage => {
            is => 'Number',
            doc => 'Desired amount of sequence in output, expressed as a multiple of the estimated genome size',
        },
    ],
    has => [
        instrument_data => {
           is => 'Genome::InstrumentData',
           is_many => 1,
           calculate => q{
              return Genome::InstrumentData->get([$self->instrument_data_id]);
           },
        },
    ],
};

sub _construct_sx_command_parts {
    my $self = shift;
    my @input_files;
    my $output_file_count = $self->get_output_file_count;
    for (my $input_number=1; $input_number <= $output_file_count; $input_number++) {
        my $new_input_file = join("/", $self->temp_staging_directory, $input_number);
        push @input_files, $new_input_file;
    }
    my $total_incoming_bases = 0;
    my @instrument_data = $self->instrument_data;
    for my $id (@instrument_data) {
        my %params = (
            instrument_data_id => $id->id,
            read_processor => $self->read_processor,
            test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
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
            return;
        }
        if ( $result->output_file_config or not defined $result->output_file_type or $result->output_file_type ne 'sanger' ) {
            $self->error_message('Incompatible output type ('.($result->output_file_config or $result->output_file_type or 'NULL').') to merge results! Can only merge sanger fastq!');
            return;
        }
        #cat the out file(s) onto the master input file(s)
        my @output_files = $result->read_processor_output_file_paths;
        for (my $input_number=0; $input_number < $output_file_count; $input_number++) {
            my $cmd = "cat ".$output_files[$input_number]." >> ".$input_files[$input_number];
            `$cmd`;
        }
        my $metrics = Genome::Model::Tools::Sx::Metrics->from_file($result->read_processor_output_metric_file);
        my $output_bases = $metrics->bases;
        $total_incoming_bases += $output_bases;
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
                return;
            }
        }

        unless (defined $estimated_genome_size) {
            $self->error_message("No estimated genome size found.  This must be specified on the instrument data if a coverage limit is needed for this result");
            return;
        }

        #convert $self->coverage and estimated_genome_size to number of bp
        my $bp_count = $self->coverage * $estimated_genome_size;
        $read_processor = "limit by-bases --incoming-sequences $total_incoming_bases --bases $bp_count";
    }

    my $output = $self->_output_config;
    unless (defined $output) {
        $self->error_message("Failed to define output");
        return;
    }
    my $input_file_type = $self->resolve_merged_input_type;
    my @inputs = map { $_.':type='.$input_file_type } @input_files;

    return ( $read_processor, \@inputs, $output );
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

sub resolve_base_name_from_instrument_data {
    my $self = shift;
    my $base = "merged";
    my @ids = $self->instrument_data_id;
    for my $id (@ids) {
        $base .= "-$id";
    }
    return $base;
}

1;

