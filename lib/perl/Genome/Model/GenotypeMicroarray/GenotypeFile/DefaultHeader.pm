package Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::GenotypeMicroarray::GenotypeFile::DefaultHeader { 
    is => 'UR::Singleton',
};

sub header_lines {
    no warnings;
    return [
        '##fileformat=VCFv4.1',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed">',
        '##INFO=<ID=CC,Number=1,Type=Float,Description="CNV Confidence">',
        '##INFO=<ID=CV,Number=1,Type=Float,Description="CNV Value">',
        '##INFO=<ID=GC,Number=1,Type=Float,Description="GC Score">',
        '##INFO=<ID=LR,Number=1,Type=Float,Description="Log R Ratio">',
        '##INFO=<ID=OG,Number=1,Type=Integer,Description="Original genotype calls">',
    ];
}

our $header;
sub header {
    return $header if $header;
    $header = Genome::File::Vcf::Header->create(lines => header_lines());
    return $header;
}

1;

