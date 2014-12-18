package Genome::VariantReporting::Suite::BamReadcount::ComponentBase;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::BamReadcountParser;

class Genome::VariantReporting::Suite::BamReadcount::ComponentBase {
    is => ['Genome::VariantReporting::Framework::Component::WithSampleName'],
    has => [
    ],
};

sub get_readcount_entries {
    my $self = shift;
    my $entry = shift;

    return Genome::File::Vcf::BamReadcountParser::get_bam_readcount_entries(
        $entry,
        $self->sample_index($entry->{header}),
    );
}

1;

