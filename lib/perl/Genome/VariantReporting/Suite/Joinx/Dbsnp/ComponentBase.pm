package Genome::VariantReporting::Suite::Joinx::Dbsnp::ComponentBase;

use strict;
use warnings;
use Genome;
use Memoize qw();

class Genome::VariantReporting::Suite::Joinx::Dbsnp::ComponentBase {
};

sub _caf_parser {
    my $self = shift;
    my $header = shift;
    return Genome::File::Vcf::DbsnpAFParser->new($header);
}

Memoize::memoize('_caf_parser', LIST_CACHE => 'MERGE');
1;

