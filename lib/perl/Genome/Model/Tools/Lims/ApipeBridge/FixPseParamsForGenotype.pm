package Genome::Model::Tools::Lims::ApipeBridge::FixPseParamsForGenotype;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Lims::ApipeBridge::FixPseParamsForGenotype { 
    is => 'Genome::Model::Tools::Lims::ApipeBridge::FixPseParamsForBase',
};

sub instrument_data_type { return 'genotyper results'; }
sub valid_prior_processes { return ( 'assign genotype qc', 'assign genotype external qc' ); }

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

