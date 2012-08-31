package Genome::Model::Event::Build::AmpliconAssembly::VerifyInstrumentData;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Event::Build::AmpliconAssembly::VerifyInstrumentData {
    is => 'Genome::Model::Event',
};

sub bsub_rusage {
    return "-R 'span[hosts=1]'";
}

sub execute {
    my $self = shift;

    if ( $self->model->sequencing_center eq 'gsc' ) {
        $self->_link_instrument_data
            or return;
    } # TODO add logic for other centers...

    my @links = glob($self->build->chromat_dir.'/*');
    
    return ( @links ) ? 1 : 0;
}

sub _link_instrument_data {
    my $self = shift;

    my @instrument_data = $self->build->instrument_data;
    unless ( @instrument_data ) {
        $self->error_message(
            sprintf(
                'No instrument data found for amplicon assembly build (<Id> %s <Model Name> %s <Model Id> %s).',
                $self->build->id,
                $self->model->name,
                $self->model->id,
            )
        );
    }

    for my $instrument_data ( @instrument_data ) {
        unless ( $instrument_data->dump_to_file_system ) {
            $self->error_message(
                sprintf(
                    'Error dumping instrument data (%s <Id> %s) for amplicon assembly build (<Id> %s <Model Name> %s <Model Id> %s).',
                    $instrument_data->run_name,
                    $instrument_data->id,
                    $self->build->id,
                    $self->model->name,
                    $self->model->id,
                )
            );
            return;
        }

        unless ( $self->build->link_instrument_data( $instrument_data ) ) {
            $self->error_message(
                sprintf(
                    'Error linking instrument data (%s <Id> %s) for amplicon assembly build (<Id> %s <Model Name> %s <Model Id> %s).',
                    $instrument_data->run_name,
                    $instrument_data->id,
                    $self->build->id,
                    $self->model->name,
                    $self->model->id,
                )
            );
            return;
        }
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
