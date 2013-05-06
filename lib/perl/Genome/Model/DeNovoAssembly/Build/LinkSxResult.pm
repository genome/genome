package Genome::Model::DeNovoAssembly::Build::LinkSxResult;

use strict;
use warnings;

use Genome;
use Genome::Model::DeNovoAssembly::SxReadProcessor;

class Genome::Model::DeNovoAssembly::Build::LinkSxResult {
    is => 'Command::V2',
    has_input => [
        build => { 
            is => 'Genome::Model::Build::DeNovoAssembly',
            is_output => 1,
            doc => 'Build to link SX result.',
        },
        sx_result => {
            is => 'Genome::InstrumentData::SxResult',
            doc => 'Sx result to link. It should have been created durin the process or merge step.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Link SX result...');

    $self->status_message('Build: '.$self->build->__display_name__);

    my $sx_result = $self->sx_result;
    $self->status_message('Sx result id: '.$sx_result->id);
    $self->status_message('Sx output dir: '.$sx_result->output_dir);
    $self->status_message('Instrument data: '.$sx_result->instrument_data_id);

    foreach my $output_file ( $sx_result->read_processor_output_files ) {
        my $target = $sx_result->output_dir.'/'.$output_file;
        my $link_name = $self->build->data_directory.'/'.$output_file;
        $self->status_message("Link output file $target to $link_name");
        Genome::Sys->create_symlink($target, $link_name);
    }

    $self->status_message('Link SX result...');
    for my $type (qw/ input output /) {
        my $metrics_file_method = 'read_processor_'.$type.'_metric_file';
        my $target = $sx_result->output_dir.'/'.$sx_result->$metrics_file_method;
        my $link_name = $self->build->data_directory.'/'.$sx_result->$metrics_file_method;
        $self->status_message("Link $type metrics file $target to $link_name");
        Genome::Sys->create_symlink($target, $link_name);
    }

    $self->status_message('Link SX result...done');
    return 1;
}

1;

