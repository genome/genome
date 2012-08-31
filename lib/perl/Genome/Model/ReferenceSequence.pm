package Genome::Model::ReferenceSequence;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence {
    is => 'Genome::ModelDeprecated',
    has => [
        prefix => {
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'prefix', value_class_name => 'UR::Value' ],
            doc => 'The source of the sequence (such as NCBI).  May not contain spaces.',
            is_mutable => 1,
            is_many => 0,
            is_optional => 1,
        },
        desc => {
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'desc', value_class_name => 'UR::Value' ],
            doc => 'Optional additional qualifier, i.e. "lite" for the GRC-human-lite refseq model.',
            is_mutable => 1,
            is_many => 0,
            is_optional => 1,
        },
        _expected_name => {
            calculate_from => ['prefix','subject_name','desc'],
            calculate => 'no warnings; my $v = ($desc ? "$prefix-$subject_name-$desc" : "$prefix-$subject_name"); $v =~ s/ /_/g; $v'
        },
        sample_names => {
            is => 'ARRAY',
            calculate => q{ return; }
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
        },
    ],
    doc => 'a versioned reference sequence, with cordinates suitable for annotation',
};

sub build_by_version {
    my $self = shift;
    my $version = shift;
    my @b = Genome::Model::Build::ImportedReferenceSequence->get(
        'version' => $version,
        'model_id' => $self->genome_model_id
    );
    if (@b > 1) {
        die "Multiple builds for version $version for model " . $self->genome_model_id;
    }
    return $b[0];
}

sub build_needed {
    my $self = shift;

    #These models always have one build per "version" of the reference.
    #The models themselves do not store information relevant to whether a future build is needed.
    return 0;
}

1;
