package Genome::ModelGroup::View::GoldSnpComparison::Xml;

use strict;
use warnings;

use Genome;

class Genome::ModelGroup::View::GoldSnpComparison::Xml {
    is => 'Genome::Model::Set::View::GoldSnpComparison::Xml',
    has_constant => [
        perspective => {
            value => 'gold-snp-comparison',
        },
    ]
};

# model groups are just sets of models, right?  lean on that
# and just provide an alternate path to getting the subjects;

sub _subject_models {
    my $self = shift;
    
    my @members = $self->subject->models;
}


1;
