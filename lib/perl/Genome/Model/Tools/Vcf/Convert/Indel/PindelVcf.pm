package Genome::Model::Tools::Vcf::Convert::Indel::PindelVcf;

use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Tools::Vcf::Convert::Indel::PindelVcf {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a TGI standard VCF file from Pindel output, TCGA-compliant or not',
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
    return 'Pindel';
}


sub get_info_meta {
    return (
        {MetaType => "INFO", ID => "END",    Number => 1,   Type => "Integer", Description => "End position of the variant described in this record"},
        {MetaType => "INFO", ID => "HOMLEN", Number => ".", Type => "Integer", Description => "Length of base pair identical micro-homology at event breakpoints"},
        {MetaType => "INFO", ID => "HOMSEQ", Number => ".", Type => "String",  Description => "Sequence of base pair identical micro-homology at event breakpoints"},
        {MetaType => "INFO", ID => "SVLEN",  Number => ".", Type => "Integer", Description => "Difference in length between REF and ALT alleles"},
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
    my @infos   = split /;/, $columns[7];
    my @new_infos;

    for my $info (@infos) {
        push @new_infos, $info unless $info =~ /^SVTYPE\=/;  #SVTYPE not useful and conflict with TCGA definition
    }

    $columns[7]   = join ';', @new_infos;
    $columns[8]  .= ':DP:BQ:SS';
    $columns[9]  .= ':.:.:.';
    $columns[10] .= ':.:.:2' if $columns[10];  #sometimes it gets germline single bam

    return join "\t", @columns;
}


1;

