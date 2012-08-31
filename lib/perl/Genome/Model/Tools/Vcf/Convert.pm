package Genome::Model::Tools::Vcf::Convert;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Vcf::Convert {
    is => ['Command'],
};

sub help_brief {
    "Tools to convert variant lists into VCF files."
}

sub help_detail {
    return <<EOS
Tools to convert variant lists into VCF files.
EOS
}

1;
