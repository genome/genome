package Genome::VariantReporting::Suite::Vep::ConsequenceFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use List::MoreUtils qw (uniq);
use Set::Scalar;

class Genome::VariantReporting::Suite::Vep::ConsequenceFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => {
        consequences => {
            is => 'Text',
            is_many => 1,
            doc => 'A list of transcript consequences',
        },
    },
    doc => q{Filter out variants whose transcript consequences don't match any of the ones in the given consequences list},
};


sub name {
    return 'consequence';
}

sub requires_annotations {
    return ('vep');
}

sub filter_entry {
    my ($self, $entry) = @_;
    my %return_values;

    my $vep_parser = Genome::File::Vcf::VepConsequenceParser->new($entry->{header});

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        my @transcripts  = $vep_parser->transcripts($entry, $alt_allele);
        my $consequences = Set::Scalar->new(uniq map {split /\&/, lc($_)} grep {defined($_)} map {$_->{consequence}} @transcripts);
        my $desired_consequences = Set::Scalar->new($self->consequences);
        if ($consequences->intersection($desired_consequences)) {
            $return_values{$alt_allele} = 1;
        }
        else {
            $return_values{$alt_allele} = 0;
        }
    }

    return %return_values;
}

sub vcf_id {
    my $self = shift;
    return 'CONSEQUENCES_' . uc(join('_', $self->consequences));
}

sub vcf_description {
    my $self = shift;
    return 'Transcript consequence is one of: ' . join(', ', $self->consequences);
}

1;




