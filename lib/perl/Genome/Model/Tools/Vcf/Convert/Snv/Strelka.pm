package Genome::Model::Tools::Vcf::Convert::Snv::Strelka;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Snv::Strelka {
    is => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from strelka output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from strelka snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}


sub source {
    return 'Strelka';
}

#maintain strelka original vcf as much as possible
sub extra_format_meta {
    return (
        {MetaType => "FORMAT", ID => "FDP",   Number => 1, Type => "Integer", Description => "Number of basecalls filtered from original read depth for tier1"},
        {MetaType => "FORMAT", ID => "SDP",   Number => 1, Type => "Integer", Description => "Number of reads with deletions spanning this site at tier1"},
        {MetaType => "FORMAT", ID => "SUBDP", Number => 1, Type => "Integer", Description => "Number of reads below tier1 mapping quality threshold aligned across this site"},
        {MetaType => "FORMAT", ID => "AU",    Number => 2, Type => "Integer", Description => "Number of 'A' alleles used in tiers 1,2"},
        {MetaType => "FORMAT", ID => "CU",    Number => 2, Type => "Integer", Description => "Number of 'C' alleles used in tiers 1,2"},
        {MetaType => "FORMAT", ID => "GU",    Number => 2, Type => "Integer", Description => "Number of 'G' alleles used in tiers 1,2"},
        {MetaType => "FORMAT", ID => "TU",    Number => 2, Type => "Integer", Description => "Number of 'T' alleles used in tiers 1,2"},
    );
}

sub get_info_meta {
    return (
        {MetaType => "INFO", ID => "QSS",     Number => 1, Type => "Integer", Description => "Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal"},
        {MetaType => "INFO", ID => "TQSS",    Number => 1, Type => "Integer", Description => "Data tier used to compute QSS"},
        {MetaType => "INFO", ID => "NT",      Number => 1, Type => "String",  Description => "Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}."},
        {MetaType => "INFO", ID => "QSS_NT",  Number => 1, Type => "Integer", Description => "Quality score reflecting the joint probability of a somatic variant and NT"},
        {MetaType => "INFO", ID => "TQSS_NT", Number => 1, Type => "Integer", Description => "Data tier used to compute QSS_NT"},
        {MetaType => "INFO", ID => "SGT",     Number => 1, Type => "String",  Description => "Most likely somatic genotype excluding normal noise states"},
    );
        ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">   leave this now
}

sub get_filter_meta {
    return (
        {MetaType => "FILTER", ID => "BCNoise", Description => "Fraction of basecalls filtered at this site in either sample is at or above 0.4"},
        {MetaType => "FILTER", ID => "SpanDel", Description => "Fraction of reads crossing site with spanning deletions in either sample exceeeds 0.75"},
        {MetaType => "FILTER", ID => "QSS_ref", Description => "Normal sample is not homozygous ref or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15"},
    );
}


sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^#/; # no vcf header here
    my @columns = split /\t/, $line;

    my ($ref, $alt, $info, $n_sample, $t_sample) = map{$columns[$_]}(3, 4, 7, 9, 10);
    my ($n_gt_info, $n_gt_str, $t_gt_str) = $info =~ /NT=(\S+?);QSS.*SGT=(\S+?)\->(\S+?);/;

    my @n_data = split /:/, $n_sample;
    my @t_data = split /:/, $t_sample;

    my @alts = split /,/, $alt;

    #sometimes ALT column gets only .
    my $n_ad = $alt eq '.' ? '.' : parse_ad(\@n_data, \@alts);
    my $t_ad = $alt eq '.' ? '.' : parse_ad(\@t_data, \@alts);

    my %ids;
    my $id = 0;

    for my $base ($ref, @alts) {
        $ids{$base} = $id;
        $id++;
    }

    my $n_gt = $n_gt_info eq 'ref' ? '0/0' : parse_gt($n_gt_str, \%ids);
    my $t_gt = parse_gt($t_gt_str, \%ids);

    $columns[7]  =~ s/SOMATIC;//;  #remove the meaningless SOMATIC, it is contained in every line
    $columns[8]  = 'GT:AD:BQ:SS:'. $columns[8];
    $columns[9]  = $n_gt . ':' . $n_ad . ':.:.:' . $n_sample;
    $columns[10] = $t_gt . ':' . $t_ad . ':.:2:' . $t_sample;

    return join "\t", @columns;
}


sub parse_ad {
    my ($data, $alts) = @_;
    my @base_cts = splice @$data, -4, 4;
    my @bases    = qw(A C G T);
    my %cts;
    my @ads;

    for my $id (0..3) {
        my ($base_ct) = $base_cts[$id] =~ /^(\d+),/;  #only take tier1 count since DP takes only tier1 too
        $cts{$bases[$id]} = $base_ct;
    }

    for my $alt (@$alts) {
        push @ads, $cts{$alt};
    }

    return join ',', @ads;
}


sub parse_gt {
    my ($gt_str, $ids) = @_;
    my @gt_ids = map{$ids->{$_}}(split //, $gt_str);
    return join '/', sort @gt_ids;
}


#EXAMPLE VCF FILE PUT OUT BY STRELKA NATIVELY
##fileformat=VCFv4.1
##fileDate=20120710
##source=strelka
##startTime=Tue Jul 10 19:35:13 2012
##content=strelka somatic snv calls
##germlineSnvTheta=0.001
##priorSomaticSnvRate=1e-06
##INFO=<ID=QSS,Number=1,Type=Integer,Description="Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal">
##INFO=<ID=TQSS,Number=1,Type=Integer,Description="Data tier used to compute QSS">
##INFO=<ID=NT,Number=1,Type=String,Description="Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.">
##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description="Quality score reflecting the joint probability of a somatic variant and NT">
##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description="Data tier used to compute QSS_NT">
##INFO=<ID=SGT,Number=1,Type=String,Description="Most likely somatic genotype excluding normal noise states">
##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth for tier1 (used+filtered)">
##FORMAT=<ID=FDP,Number=1,Type=Integer,Description="Number of basecalls filtered from original read depth for tier1">
##FORMAT=<ID=SDP,Number=1,Type=Integer,Description="Number of reads with deletions spanning this site at tier1">
##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description="Number of reads below tier1 mapping quality threshold aligned across this site">
##FORMAT=<ID=AU,Number=2,Type=Integer,Description="Number of 'A' alleles used in tiers 1,2">
##FORMAT=<ID=CU,Number=2,Type=Integer,Description="Number of 'C' alleles used in tiers 1,2">
##FORMAT=<ID=GU,Number=2,Type=Integer,Description="Number of 'G' alleles used in tiers 1,2">
##FORMAT=<ID=TU,Number=2,Type=Integer,Description="Number of 'T' alleles used in tiers 1,2">
##FILTER=<ID=DP,Description="Greater than 3.0x chromosomal mean depth in Normal sample">
##FILTER=<ID=BCNoise,Description="Fraction of basecalls filtered at this site in either sample is at or above 0.4">
##FILTER=<ID=SpanDel,Description="Fraction of reads crossing site with spanning deletions in either sample exceeeds 0.75">
##FILTER=<ID=QSS_ref,Description="Normal sample is not homozygous ref or ssnv Q-score < 15, ie calls with NT!=ref or QSS_NT < 15">
##cmdline=/gscmnt/gc2142/techd/tools/strelka/v0.4.6.2/strelka_workflow/strelka/scripts/consolidateResults.pl --config=/gscmnt/gc2142/techd/analysis/strelka/all_chrs/results/config/run.config.ini
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
#1       10231   .       C       A       .       QSS_ref NT=ref;QSS=1;QSS_NT=1;SGT=AC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    32:4:8:0:0,3:28,60:0,0:0,1      84:6:69:0:7,21:71,192:0,0:0,1
#1       10333   .       C       T       .       QSS_ref NT=ref;QSS=4;QSS_NT=4;SGT=CC->CT;SOMATIC;TQSS=1;TQSS_NT=1       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    13:1:0:0:0,0:12,33:0,0:0,1      49:8:2:0:0,0:37,92:0,0:4,10
#1       10440   .       C       A       .       QSS_ref NT=ref;QSS=5;QSS_NT=5;SGT=CC->AC;SOMATIC;TQSS=2;TQSS_NT=2       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    16:1:27:0:0,1:15,56:0,0:0,0     107:2:69:0:11,18:94,204:0,0:0,0
#1       10473   .       G       A       .       QSS_ref NT=ref;QSS=6;QSS_NT=6;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1       DP:FDP:SDP:SUBDP:AU:CU:GU:TU    9:1:0:0:0,0:0,1:8,13:0,0        62:1:0:0:18,23:0,2:43,52:0,0



