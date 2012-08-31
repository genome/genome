package Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForSolexa;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForSolexa { 
    is => 'Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForBase',
};

sub instrument_data_type { return 'solexa'; }
sub valid_prior_processes { return ( 'copy sequence files' ); }

sub _get_sequence_item_from_prior {
    my ($self, $prior) = @_;

    my $index_illumina = $prior->get_index_illumina;
    if ( not $index_illumina ) {
        $self->error_message('No index illumin for prior PSE!');
        next;
    }
    $self->status_message('Index illumina: '.$index_illumina->id);

    return$index_illumina;
}

1;

