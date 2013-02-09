package Genome::Model::ReferenceSequence;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceSequence {
    is => 'Genome::ModelDeprecated',
    has => [
        source => {
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'prefix', value_class_name => 'UR::Value' ],
            doc => 'The source of the sequence (such as GRC).  May not contain spaces.',
            is_mutable => 1,
            is_many => 0,
            is_optional => 1,
        },
        prefix => {
            is_deprecated => 1,
            is => 'UR::Value',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'prefix', value_class_name => 'UR::Value' ],
            doc => 'The source of the sequence (such as GRC).  May not contain spaces.',
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
    if ($version eq '36' or $version eq '37') {
        # this is present only to help developers troubleshoot problems related to poorly maintained refseq data
        # remove it when this data is cleaned-up (ssmith-2013-02-09)
        Carp::cluck("Cearching for version $version is likely due to out-of-date code.  Reference versions now end in -lite, or a patch level indicator in most cases");
    }
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
