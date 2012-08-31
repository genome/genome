package Genome::Model::SomaticCapture::Report::Variant;

use strict;
use warnings;

use Genome;

#For now, just use the Somatic pipeline's Variant report.
class Genome::Model::SomaticCapture::Report::Variant {
    is => 'Genome::Model::Somatic::Report::Variant',
    doc => 'A summary of predicted somatic variants for the model.',
};

1;
