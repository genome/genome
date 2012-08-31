package Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForGenotype;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForGenotype { 
    is => 'Genome::Model::Tools::Lims::ApipeBridge::FixPidfaParamsForBase',
    has_optional => [
        instrument_data_id => {
            is => 'Integer',
            doc => 'Instrument data id to get PIDFA to fix',
        },
    ],
};

sub instrument_data_type { return 'genotyper results'; }
sub valid_prior_processes { return ( 'assign genotype qc', 'assign genotype external qc' ); }
sub _additional_starting_points { return 'instrument_data_id'; }

sub _init_with_instrument_data_id {
    my $self = shift;

    my @pse_params = GSC::PSEParam->get(
        param_name => 'instrument_data_id',
        param_value => $self->instrument_data_id,
    );
    if ( not @pse_params ) { 
        $self->error_message('Failed to get PSE params for instrument data id! '.$self->instrument_data_id);
        return;
    }

    my @pidfas = GSC::PSE->get(id => [ map { $_->pse_id } @pse_params ], ps_id => 3870);
    if ( not @pidfas ) {
        $self->error_message('Failed to get PIDFA for pse param ids! '.join(' ', map { $_->id } @pidfas));
        return;
    }
    elsif ( @pidfas > 1 ) {
        $self->error_message('Got '.@pidfas.' PIDFAs for instrument data id! '.$self->pidfa_id);
        return;
    }

    my ($tp_pse) = GSC::TppPSE->get(pse_id => $pidfas[0]->id);
    if ( not $tp_pse ) {
        $self->error_message('No transfer pattern pse for PIDFA! '.$pidfas[0]->id);
        return;
    }#my ($prior_id) = $pidfa->added_param('control_pse_id'); 

    if ( not $tp_pse->prior_pse_id ) {
        $self->error_message('No prior pse id iun transfer pattern for PIDFA!');
        return;
    }

    my $prior = GSC::PSE->get(pse_id => $tp_pse->prior_pse_id);
    if ( not $prior ) {
        $self->error_message('Failed to get prior PSE for id!'. $tp_pse->prior_pse_id);
        return;
    }

    return ($pidfas[0], $prior);
}

sub _get_sequence_item_from_prior {
    my ($self, $prior) = @_;

    my @genotypes = eval{ $prior->get_genotype; };
    @genotypes = eval{ $prior->get_external_genotype; } if not @genotypes;
    if ( not @genotypes ) {
        $self->error_message("No genotypes for prior PSE! ".$prior->id);
    }
    elsif ( @genotypes > 1 ) {
        $self->error_message("Multiple genotypes for prior PSE! ".$prior->id);
    }

    if ( $genotypes[0]->status and $genotypes[0]->status ne 'pass' ) {
        $self->error_message("Invalid status for genotype! ".$genotypes[0]->id.' '.$genotypes[0]->status);
        return;
    }

    return $genotypes[0];
}

sub _additional_params_to_fix {
    my ($self, $params_to_fix, $sequence_item) = @_;
    $params_to_fix->{genotype_file} = $sequence_item->get_genotype_file;
    return 1;
}

1;

