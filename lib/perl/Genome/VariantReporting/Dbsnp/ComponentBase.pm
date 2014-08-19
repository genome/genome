package Genome::VariantReporting::Dbsnp::ComponentBase;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Dbsnp::ComponentBase {
};

sub _caf_parser {
    my $self = shift;
    my $header = shift;
    return Genome::File::Vcf::DbsnpAFParser->new($header);
}

Memoize::memoize('_caf_parser');
1;

