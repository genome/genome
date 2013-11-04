package Genome::Model::Tools::Vcf::Convert::Indel::VarscanSomaticValidation;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use POSIX qw(log10);

class Genome::Model::Tools::Vcf::Convert::Indel::VarscanSomaticValidation {
    is => 'Genome::Model::Tools::Vcf::Convert::Indel::VarscanSomatic',
    doc => 'Generate a VCF file from varscan somatic validation output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from varscan somatic validation indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    my $self = shift;
    return "VarscanSomaticValidation";
}
