package Genome::Site::TGI::SnpConcordance;

use strict;
use warnings;
use Genome;

class Genome::Site::TGI::SnpConcordance {
    table_name => 'SNP_CONCORDANCE',
    data_source => 'Genome::DataSource::Dwrac',
    id_by => [
        creation_event_id => { is => 'Number' },
    ],
    has => [
        seq_id => { is => 'Number' },
        replicate_seq_id => { is => 'Number' },
    ],
    has_optional => [
        no_data => { is => 'Number' },
        match => { is => 'Number' },
        no_data_1 => { is => 'Number' },
        no_data_2 => { is => 'Number' },
        reverse_complemented => { is => 'Number' },
        loh => { is => 'Number' },
        not_match => { is => 'Number' },
        total => { is => 'Number' },
    ],
};


sub match_percent {
    my $self = shift;
    
    my $match = $self->match;
    my $total = $self->total;

    my $match_percent;
    if (defined $match && $total) { # match could be zero but total shouldn't
        $match_percent = $match / $total * 100;
    }

    return $match_percent;
}


sub is_internal_comparison {
    my $self = shift;
    my $replicate_seq_id = $self->replicate_seq_id;

    my $genotype;
    if ($replicate_seq_id) {
        $genotype = Genome::Site::TGI::IlluminaGenotyping->get(seq_id => $replicate_seq_id);
    }

    return defined $genotype;
}


sub is_external_comparison {
    my $self = shift;
    my $replicate_seq_id = $self->replicate_seq_id;

    my $genotype;
    if ($replicate_seq_id) {
        $genotype = Genome::Site::TGI::ExternalGenotyping->get(seq_id => $replicate_seq_id);
    }

    return defined $genotype;
}


1;

