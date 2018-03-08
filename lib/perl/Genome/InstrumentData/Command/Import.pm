package Genome::InstrumentData::Command::Import;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import {
    is => 'Command',
    doc => 'import external instrument data'
};

sub validate_target_region {
    my $self = shift;
    my $target_region = $self->target_region;
    return 0 unless defined $target_region;

    my @feature_lists = Genome::FeatureList->get(name => $target_region);
    if (not @feature_lists or @feature_lists > 1) {
        return 0;
    }
    return 1;
}

1;
