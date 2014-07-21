package Genome::VariantReporting::Reporter::WithHeaderAndSampleNames;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Reporter::WithHeaderAndSampleNames {
    is => [ 'Genome::VariantReporting::Reporter::WithHeader', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    is_abstract => 1,
    has => [
    ],
};


sub available_fields_for_interpreter {
    my $self = shift;
    my $interpreter = shift;

    return $interpreter->available_fields($self->sample_names);
}

1;
