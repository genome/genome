package Genome::Model::DeNovoAssembly::Build::MergeAndLinkSxResults;

use strict;
use warnings;

use Genome;

use Genome::Model::DeNovoAssembly::SxReadProcessor;
use Data::Dumper 'Dumper';

class Genome::Model::DeNovoAssembly::Build::MergeAndLinkSxResults {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_many => 1,
            doc => 'Build to process instrument data.',
        },
    ],
    has_output => [
        output_build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            doc => 'Build to process instrument data.',
        },
        sx_results => {
            is => 'Genome::InstrumentData::SxResult',
            is_many => 1,
            doc => 'Sx results from process instrument data.',
        },
    ],
    has_constant => [
        lsf_resource => { 
            is => 'Text',
            default_value => "-R 'select[type==LINUX64 && gtmp>1000] rusage[gtmp=1000] span[hosts=1]'",
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Merge and Link SX results...');

    my ($build) = $self->build;
    $self->status_message('Build: '.$build->__display_name__);
    $build->reads_attempted(0);
    $build->reads_processed(0);

    my $read_processor = $build->processing_profile->read_processor;
    $self->status_message('Read processor: '.$read_processor);

    my $sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
        processor => $read_processor, # what if no processing is desired?
    );
    if ( not $sx_processor ) {
        $self->error_message('Failed to create SX read processor!');
        return;
    }
    $sx_processor->determine_processing( $build->instrument_data );

    my @final_sx_result_params = $sx_processor->final_sx_result_params;
    return if not @final_sx_result_params;

    my @sx_results;
    for my $sx_result_params ( @final_sx_result_params ) {
        $self->status_message('Sx result params: '.Dumper($sx_result_params));
        my $sx_result;
        if ( $sx_result_params->{coverage} ) { # merged: get, then get_or_create
            $sx_result = Genome::InstrumentData::MergedSxResult->get_with_lock(%$sx_result_params);
            if ( not $sx_result ) {
                $sx_result = Genome::InstrumentData::MergedSxResult->get_or_create(%$sx_result_params);
            }
        }
        else { # regular: get only
            $sx_result = Genome::InstrumentData::SxResult->get_with_lock(%$sx_result_params);
        }
        return if not $sx_result;

        $self->status_message('Sx result id: '.$sx_result->id);
        $self->status_message('Sx output dir: '.$sx_result->output_dir);
        $self->status_message('Instrument data: '.join(' ', $sx_result->instrument_data_id));

        my $link_ok = $self->_link_sx_result($sx_result);
        return if not $link_ok;

        my $collect_metrics = $self->_collect_metrics_from_sx_result($sx_result);
        return if not $collect_metrics;

        push @sx_results, $sx_result;
    }

    $self->status_message('Reads attempted: '.$build->reads_attempted);
    $self->status_message('Reads processed: '.$build->reads_processed);
    $build->reads_processed_success( $build->reads_attempted ? sprintf('%0.3f', $build->reads_processed / $build->reads_attempted) : 0);
    $self->status_message('Reads processed success: '.$build->reads_processed_success);

    $self->output_build($build);
    $self->sx_results(\@sx_results);
    $self->status_message('Merge and Link SX results...OK');
    return 1;
}

sub _link_sx_result {
    my ($self, $sx_result) = @_;
    $self->status_message('Link SX result...');

    Carp::confess('No SX result given to link!') if not $sx_result;

    my ($build) = $self->build;
    foreach my $output_file ( $sx_result->read_processor_output_files ) {
        my $target = $sx_result->output_dir.'/'.$output_file;
        my $link_name = $build->data_directory.'/'.$output_file;
        $self->status_message("Link output file $target to $link_name");
        Genome::Sys->create_symlink($target, $link_name);
    }

    $self->status_message('Link SX result...');
    for my $type (qw/ input output /) {
        my $metrics_file_method = 'read_processor_'.$type.'_metric_file';
        my $target = $sx_result->$metrics_file_method;
        my $metrics_file_base_name_method = 'read_processor_'.$type.'_metric_file_base_name';
        my $link_name = $build->data_directory.'/'.$sx_result->$metrics_file_base_name_method;
        $self->status_message("Link $type metrics file $target to $link_name");
        Genome::Sys->create_symlink($target, $link_name);
    }

    $self->status_message('Link SX result...done');
    return 1;
}

sub _collect_metrics_from_sx_result {
    my ($self, $sx_result) = @_;
    $self->status_message('Collect metrics from SX result...');

    Carp::confess('No SX result given to collect metrics!') if not $sx_result;

    my ($build) = $self->build;
    my %types_and_values = ( input => 0, output => 0, );
    for my $type ( keys %types_and_values ) {
        my $metrics_file_method = 'read_processor_'.$type.'_metric_file';
        my $metrics_file = $sx_result->$metrics_file_method;
        if ( not $metrics_file or not -s $metrics_file ) {
            $self->error_message('No output metrics file for SX result! '.$sx_result->id);
            return;
        }
        $self->status_message(ucfirst($type).' metrics file: '.$metrics_file);

        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->from_file($metrics_file);
        if ( not $metrics ) {
            $self->error_message('Failed to create SX metrics from file! '.$metrics_file);
            return;
        }

        my $count = $metrics->count;
        if ( not defined $count ) {
            $self->error_message('No count in metrics! '.Data::Dumper::Dumper($metrics));
            return;
        }

        $self->status_message('Count: '.$count);
        $types_and_values{$type} = $count;
    }

    $build->reads_attempted( $build->reads_attempted + $types_and_values{input} );
    $build->reads_processed( $build->reads_processed + $types_and_values{output} );

    $self->status_message('Collect metrics from SX result...OK');
    return 1;
}

1;

