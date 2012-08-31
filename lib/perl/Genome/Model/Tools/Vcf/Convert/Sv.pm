package Genome::Model::Tools::Vcf::Convert::Sv;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::Convert::Sv {
    is => ['Command'],
};

sub help_brief {
    "Tools to convert lists of svs into VCF files."
}

sub help_detail {
    return <<EOS
Tools to convert lists of svs into VCF files.
EOS
}

1;
