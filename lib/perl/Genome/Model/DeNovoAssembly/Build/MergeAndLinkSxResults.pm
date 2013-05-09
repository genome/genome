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
            calculate_from => [qw/ input_build /],
            calculate => sub { return $_[0]; },
        },
        sx_results => {
            is => 'Genome::InstrumentData::SxResult',
            is_many => 1,
        },
    ],
    has_constant => [
        lsf_resource => { 
            is => 'Text',
            default_value => "-R 'select[type==LINUX64 && tmp>200000] rusage[tmp=200000] span[hosts=1]'",
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Merge and Link SX results...');

    my ($build) = $self->build;
    $self->status_message('Build: '.$build->__display_name__);
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

        my $link_ok = $self->_link_sx_result($sx_result);
        return if not $link_ok;

        push @sx_results, $sx_result;
    }

    $self->sx_results(\@sx_results);
    $self->status_message('Merge and Link SX results...');
    return 1;
}

sub _link_sx_result {
    my ($self, $sx_result) = @_;
    $self->status_message('Link SX result...');

    $self->status_message('Sx result id: '.$sx_result->id);
    $self->status_message('Sx output dir: '.$sx_result->output_dir);
    $self->status_message('Instrument data: '.$sx_result->instrument_data_id);

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
        my $target = $sx_result->output_dir.'/'.$sx_result->$metrics_file_method;
        my $link_name = $build->data_directory.'/'.$sx_result->$metrics_file_method;
        $self->status_message("Link $type metrics file $target to $link_name");
        Genome::Sys->create_symlink($target, $link_name);
    }

    $self->status_message('Link SX result...done');
    return 1;
}

1;

