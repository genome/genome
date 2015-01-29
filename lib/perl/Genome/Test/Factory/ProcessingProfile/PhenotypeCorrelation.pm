package Genome::Test::Factory::ProcessingProfile::PhenotypeCorrelation;

use base Genome::Test::Factory::ProcessingProfile;

our @required_params = qw(alignment_strategy trait_type cohort_type);

sub create_alignment_strategy {
    return 'instrument_data aligned to reference_sequence_build using bwa 0.5.9';
}

sub create_trait_type {
    return 'case-control';
}

sub create_cohort_type {
    return 'unrelated';
}

1;
