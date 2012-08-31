package Genome::Model::Tools::Vcf::Convert::Maf;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::Convert::Maf {
    is => ['Command'],
};

sub help_brief {
    "Tools to convert VCF files into Maf files."
}

sub help_detail {
    return <<EOS
Tools to convert VCF files into Maf files.
EOS
}

1;
