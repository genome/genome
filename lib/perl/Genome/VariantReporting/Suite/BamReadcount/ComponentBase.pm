package Genome::VariantReporting::Suite::BamReadcount::ComponentBase;

use strict;
use warnings;
use Genome;
use Genome::File::BamReadcount::Entry;
use Genome::File::Vcf::BamReadcountParser;

class Genome::VariantReporting::Suite::BamReadcount::ComponentBase {
    is => ['Genome::VariantReporting::Framework::Component::WithSampleName'],
    has => [
    ],
};

sub get_readcount_entry {
    my $self = shift;
    my $entry = shift;

    my $bam_readcount_string = $entry->sample_field($self->sample_index($entry->{header}), 'BRCT');
    return unless $bam_readcount_string and $bam_readcount_string ne '.';
    return Genome::File::BamReadcount::Entry->new(
        Genome::File::Vcf::BamReadcountParser::decode($bam_readcount_string));
}

1;

