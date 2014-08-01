package Genome::VariantReporting::Vep::MaxSegdupFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;

class Genome::VariantReporting::Vep::MaxSegdupFilter {
    is  => ['Genome::VariantReporting::Framework::Component::Filter'],
    has => [
        info_tag => {
            is  => 'String',
            doc => 'custom tag name in the info field to show segment duplication percentage, like SEGDUP=82',
        },
        threshold => {
            is  => 'Number',
            doc => 'The maximum segdup thresold percentage to pass',
        }
    ],
};

sub name {
    return 'max-segdup';
}

sub requires_annotations {
    return ('vep');
}

sub filter_entry {
    my ($self, $entry) = @_;
    my $threshold = $self->threshold;

    unless (defined $threshold and $threshold =~ /^\d+\.?\d*$/) {
        die "the input threshold: $threshold is not a number\n";
    }

    unless (exists $entry->{header}->info_types->{$self->info_tag}) {
        die sprintf("VCF header does not contain info tag (%s)", $self->info_tag);
    }

    my %return_values;

    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        my $segdup_percent = $entry->info($self->info_tag);
        $segdup_percent =~ s/\%$//;
        $return_values{$alt_allele} = $segdup_percent > $self->threshold ? 0 : 1;
    }

    return %return_values;
}


1;
