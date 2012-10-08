package Genome::Model::Tools::Vcf::Convert::Indel::PindelTcga;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::PindelTcga {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a TCGA-compliant VCF file from Pindel output',
    has => [
        column_header => {
            type => 'String',
            doc  => 'Provide vcf column header',
        },
    ],
};


sub help_synopsis {
    <<'HELP';
    Generate a TCGA-compliant VCF file from Pindel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a TCGA-complinat VCF containing all the required fields.
HELP
}

sub source {
    return 'PindelTcga';
}


sub _get_header_columns {
    return split /\t/, shift->column_header;
}


sub get_info_meta {
    return (
        {MetaType => "INFO", ID => "END",    Number => 1,   Type => "Integer", Description => "End position of the variant described in this record"},
        {MetaType => "INFO", ID => "HOMLEN", Number => 1,   Type => "Integer", Description => "Length of base pair identical micro-homology at event breakpoints"},
        {MetaType => "INFO", ID => "HOMSEQ", Number => ".", Type => "String",  Description => "Sequence of base pair identical micro-homology at event breakpoints"},
        {MetaType => "INFO", ID => "SVLEN",  Number => 1,   Type => "Integer", Description => "Difference in length between REF and ALT alleles"},
        {MetaType => "INFO", ID => "SVTYPE", Number => 1,   Type => "String",  Description => "Type of structural variant"},
        {MetaType => "INFO", ID => "NTLEN",  Number => ".", Type => "Integer", Description => "Number of bases inserted in place of deleted code"},
    );
}


#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	normal	tumor
#1	111635011	.	C	CCT	.	PASS	END=111635012;HOMLEN=3;HOMSEQ=CTC;SVLEN=2;SVTYPE=INS	GT:AD	0/1:1	0/1:3
#1	121013459	.	A	AT	.	PASS	END=121013460;HOMLEN=3;HOMSEQ=TTT;SVLEN=1;SVTYPE=INS	GT:AD	0/1:12	0/1:5
#1	178994538	.	G	GGCCCCTGCGCA	.	PASS	END=178994539;HOMLEN=0;SVLEN=11;SVTYPE=INS	GT:AD	0/1:3	0/1:4

sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^#/;  #skip the header
    my @columns = split /\t/, $line;

    $columns[8]  .= ':DP:BQ:SS';
    $columns[9]  .= ':.:.:.';
    $columns[10] .= ':.:.:2';

    return join "\t", @columns;
}


1;

