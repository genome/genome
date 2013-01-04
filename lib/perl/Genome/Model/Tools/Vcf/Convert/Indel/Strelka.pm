package Genome::Model::Tools::Vcf::Convert::Indel::Strelka;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Indel::Strelka {
    is => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from strelka output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from strelka indel output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the indels.
HELP
}

sub source {
    return "Strelka";
}

sub common_format_meta {
    my $self = shift;
    return (
        $self->SUPER::common_format_meta,
        #using the common FORMAT DP tag description to avoid merge conflict. strelka DP is for tier1 only
        {MetaType => "FORMAT", ID => "DP2",     Number => 1, Type => "Integer", Description => "Read depth for tier2"},
        {MetaType => "FORMAT", ID => "TAR",     Number => 2, Type => "Integer", Description => "Reads strongly supporting alternate allele for tiers 1,2"},
        {MetaType => "FORMAT", ID => "TIR",     Number => 2, Type => "Integer", Description => "Reads strongly supporting indel allele for tiers 1,2"},
        {MetaType => "FORMAT", ID => "TOR",     Number => 2, Type => "Integer", Description => "Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2"},
        {MetaType => "FORMAT", ID => "DP50",    Number => 1, Type => "Float",   Description => "Average tier1 read depth within 50 bases"},
        {MetaType => "FORMAT", ID => "FDP50",   Number => 1, Type => "Float",   Description => "Average tier1 number of basecalls filtered from original read depth within 50 bases"},
        {MetaType => "FORMAT", ID => "SUBDP50", Number => 1, Type => "Float",   Description => "Average number of reads below tier1 mapping quality threshold aligned across sites within 50 bases"},
        {MetaType => "FORMAT", ID => "FT",      Number => 1, Type => "String",  Description => "Sample genotype filter"},
    );
}


sub get_info_meta {
    return (
        {MetaType => "INFO", ID => "QSI",     Number => 1, Type => "Integer", Description => "Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal"},
        {MetaType => "INFO", ID => "TQSI",    Number => 1, Type => "Integer", Description => "Data tier used to compute QSI"},
        {MetaType => "INFO", ID => "NT",      Number => 1, Type => "String",  Description => "Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}."},
        {MetaType => "INFO", ID => "QSI_NT",  Number => 1, Type => "Integer", Description => "Quality score reflecting the joint probability of a somatic variant and NT"},
        {MetaType => "INFO", ID => "TQSI_NT", Number => 1, Type => "Integer", Description => "Data tier used to compute QSI_NT"},
        {MetaType => "INFO", ID => "SGT",     Number => 1, Type => "String",  Description => "Most likely somatic genotype excluding normal noise states"},
        {MetaType => "INFO", ID => "RU",      Number => 1, Type => "String",  Description => "Smallest repeating sequence unit in inserted or deleted sequence"},
        {MetaType => "INFO", ID => "RC",      Number => 1, Type => "Integer", Description => "Number of times RU repeats in the reference allele"},
        {MetaType => "INFO", ID => "IC",      Number => 1, Type => "Integer", Description => "Number of times RU repeats in the indel allele"},
        {MetaType => "INFO", ID => "IHP",     Number => 1, Type => "Integer", Description => "Largest reference interupted homopolymer length intersecting with the indel"},
        {MetaType => "INFO", ID => "SVTYPE",  Number => 1, Type => "String",  Description => "Type of structural variant"},
        {MetaType => "INFO", ID => "OVERLAP", Number => 0, Type => "Flag",    Description => "Somatic indel possibly overlaps a second indel."},
    );
}


sub get_filter_meta {
    return (
        {MetaType => "FILTER", ID => "DP",      Description => "Greater than 3.0x chromosomal mean depth in Normal sample"},
        {MetaType => "FILTER", ID => "Repeat",  Description => "Sequence repeat of more than 8x in the reference sequence"},
        {MetaType => "FILTER", ID => "iHpol",   Description => "Indel overlaps an interupted homopolymer longer than 14x in the reference sequence"},
        {MetaType => "FILTER", ID => "BCNoise", Description => "Average fraction of filtered basecalls within 50 bases of the indel exceeds 0.3"},
        {MetaType => "FILTER", ID => "QSI_ref", Description => "Normal sample is not homozygous ref or sindel Q-score < 30, ie calls with NT!=ref or QSI_NT < 30"},
    );
}


sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^#/; # no vcf header here
    my @columns = split /\t/, $line;

    my ($filter, $info, $n_sample, $t_sample) = map{$columns[$_]}(6, 7, 9, 10);
    my ($n_gt_info, $t_gt_info) = $info =~ /;NT=(\S+?);.*SGT.*\->(\S+?);/;

    my %gt_info = (
        ref => '0/0',
        het => '0/1',
        hom => '1/1',
        conflict => '.',
    );

    my $n_gt = $gt_info{$n_gt_info};
    my $t_gt = $gt_info{$t_gt_info};

    my ($n_ad) = $n_sample =~ /^\d+:\d+:(\d+),/;
    my ($t_ad) = $t_sample =~ /^\d+:\d+:(\d+),/;

    $columns[7]  =~ s/SOMATIC;//;  #remove the meaningless SOMATIC, it is contained in every line
    $columns[8]  = 'GT:AD:BQ:SS:'. $columns[8] . ':FT';
    $columns[9]  = $n_gt . ':' . $n_ad . ':.:.:' . $n_sample . ':' . $filter;
    $columns[10] = $t_gt . ':' . $t_ad . ':.:2:' . $t_sample . ':' . $filter;

    return join "\t", @columns;
}


#sample original strelka indel vcf output
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	H_JG-701414_763290-1205725	H_JG-701414_763290-1206178
#1	965051	.	ATGTGTG	A	.	QSI_ref	IC=5;IHP=2;NT=ref;QSI=1;QSI_NT=1;RC=8;RU=TG;SGT=het->het;SOMATIC;TQSI=1;TQSI_NT=1	DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50	8:8:6,6:0,0:2,4:10.3:0.00:0.00	18:18:8,8:5,6:5,8:21:0.25:0.00
#1	1034700	.	TTGGAGA	T	.	QSI_ref	IC=0;IHP=2;NT=ref;QSI=16;QSI_NT=16;RC=1;RU=TGGAGA;SGT=ref->het;SOMATIC;TQSI=1;TQSI_NT=1	DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50	39:39:39,49:0,0:2,4:47.77:1.34:0.00	74:74:64,79:3,4:10,17:83.09:2.64:0.00
#1	2212562	.	G	GA	.	QSI_ref	IC=1;IHP=3;NT=ref;QSI=10;QSI_NT=10;RC=0;RU=A;SGT=ref->het;SOMATIC;TQSI=1;TQSI_NT=1	DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50	12:12:12,18:0,0:0,0:12.17:0.00:0.00	14:14:10,22:3,3:1,3:13.93:0.00:0.00


