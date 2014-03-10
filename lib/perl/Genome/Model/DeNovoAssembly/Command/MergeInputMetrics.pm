package Genome::Model::DeNovoAssembly::Command::MergeInputMetrics;

use strict;
use warnings;

use Genome;

require File::Temp;

class Genome::Model::DeNovoAssembly::Command::MergeInputMetrics {
    is => 'Command::V2',
     has_optional => [
        _input_count => { is => 'Number', default_value => 0, },
        _input_bases => { is => 'Number', default_value => 0, },
        _output_count => { is => 'Number', default_value => 0, },
        _output_bases => { is => 'Number', default_value => 0, },
        _original_base_limit => { is => 'Number', },
        _base_limit => { is => 'Number', },
    ],

    has_input => [
        input_builds => { is => 'Genome::Model::Build::DeNovoAssembly',
            is_many => 1},
    ],
    has_output => [
        output_build => { is => 'Genome::Model::Build::DeNovoAssembly',
            calculate_from => ['input_builds'],
            calculate => sub { return $_[0]; }
        },
    ],
};

sub execute {
    my $self = shift;

    my @instrument_data = $self->output_build->instrument_data;
    #TODO: add in base_limit stuff here as well
    #$self->_setup_base_limit;

    $self->debug_message('Merge instrument data');
    INST_DATA: for my $instrument_data (@instrument_data) {
        $self->_update_metrics($instrument_data);
            
        #last INST_DATA if $self->_has_base_limit_been_reached;
    }

    $self->debug_message('Merge instrument data...OK');

    my $reads_attempted = $self->_input_count;
    my $reads_processed = $self->_output_count;
    my $reads_processed_success = ( $reads_attempted ? sprintf('%0.3f', $reads_processed / $reads_attempted) : 0);
    $self->output_build->add_metric(name => 'reads_attempted', value => $reads_attempted);
    $self->output_build->add_metric(name => 'reads_processed', value => $reads_processed);
    $self->output_build->add_metric(name => 'reads_processed_success', value => $reads_processed_success);
    $self->debug_message('Reads attempted: '.$reads_attempted);
    $self->debug_message('Reads processed: '.$reads_processed);
    $self->debug_message('Reads processed success: '.($reads_processed_success * 100).'%');

    $self->debug_message('Merge instrument data...OK');
    return 1;
}

sub _setup_base_limit {
    my $self = shift;

    my $base_limit = $self->output_build->calculate_base_limit_from_coverage;
    return 1 if not defined $base_limit;

    $self->debug_message('Setting base limit to: '.$base_limit);
    $self->_original_base_limit($base_limit);
    $self->_base_limit($base_limit);

    return 1;
}

sub _update_metrics {
    my $self = shift;
    my $instrument_data = shift;
    $self->debug_message('Update metrics...');

    for my $type (qw/ input output /) {
        my $metrics_file_method = $type.'_metrics_file_for_instrument_data';
        my $metrics_file = $self->output_build->$metrics_file_method($instrument_data);
        $self->debug_message(ucfirst($type)." file: $metrics_file");
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
            $self->debug_message("Update $type $name from $metric to $new_metric");
        }
    }

    $self->debug_message('Update metrics...OK');
    return 1;
}

sub _has_base_limit_been_reached {
    my $self = shift;

    return if not defined $self->_base_limit;

    $self->debug_message('Original base limit: '.$self->_original_base_limit);
    $self->debug_message('Bases processed: '.$self->_output_bases);
    my $current_base_limit = $self->_original_base_limit - $self->_output_bases;
    $self->_base_limit($current_base_limit);
    if ( $current_base_limit <= 0 ) {
        $self->debug_message('Reached base limit. Stop processing!');
        return 1;
    }
    $self->debug_message('New base limit: '.$self->_base_limit);

    $self->debug_message('Base limit not reached. Continue processing.');
    return;
}

1;
