package Genome::Model::Tools::Vcf::Convert::Snv;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::Convert::Snv {
    is => ['Command'],
};

sub help_brief {
    "Tools to convert lists of snvs into VCF files."
}

sub help_detail {
    return <<EOS
Tools to convert lists of snvs into VCF files.
EOS
}

1;
