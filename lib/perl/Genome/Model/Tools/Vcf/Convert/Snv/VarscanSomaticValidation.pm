package Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;
use POSIX qw(log10);

class Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomaticValidation {
    is => 'Genome::Model::Tools::Vcf::Convert::Snv::VarscanSomatic',
    doc => 'Generate a VCF file from varscan somatic validationoutput'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from varscan somatic validation snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "VarscanSomaticValidation";
}
