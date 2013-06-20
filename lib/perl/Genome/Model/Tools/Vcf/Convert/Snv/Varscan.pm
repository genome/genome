package Genome::Model::Tools::Vcf::Convert::Snv::Varscan;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Snv::Varscan {
    is =>  'Genome::Model::Tools::Vcf::Convert::Base' ,
    doc => 'Generate a VCF file from varscan output'
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from varscan snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    my $self = shift;
    return "Varscan";
}

sub parse_line { 
    my $self=shift;
    my $line = shift;

    my @fields = split "\t", $line;
    my $chr = $fields[0];
    my $pos = $fields[1];

    my $ref_allele = $fields[2];
    $ref_allele =~ s/[^ACGTN\-]/N/g;
    my $var_allele_iub = $fields[3];
    #convert iub to only variant alleles
    my @var_alleles =   Genome::Info::IUB->variant_alleles_for_iub($ref_allele, $var_allele_iub);
    my @all_alleles = Genome::Info::IUB->iub_to_alleles($var_allele_iub);
    my $GT = $self->generate_gt($ref_allele, \@var_alleles, \@all_alleles);
    if(!$GT) {
        $self->error_message("unable to convert $line into GT field");
        return 0;
    }
    my $alt_alleles = join(",", @var_alleles);
    # We infer a Phred-like quality from the VarScan P-value
    my $GQ= 0 - (10 * log($fields[11]) / log(10)) if($fields[11] && $fields[11] != 0);
    $GQ = sprintf("%d", $GQ);
    $GQ = 255 if($GQ > 255);
    #total depth = "reads1 + reads2"
    my $DP= $fields[4]+$fields[5];
    #avg base quality for var = avg_qual2. For Ref it would be $fields[9] which is avg_qual1
    my $refBQ = $fields[9];
    my $BQ = join(",",($refBQ, $fields[10]));
    #avg mapping quality. Previously this next line read my MQ = $fields[13]; This is not output by VarScan so we should hardcode to '.'
    my $MQ = ".";
    #allele_depth = "reads2"
    my $refAD = $fields[4];
    my $AD = join(",",($refAD, $fields[5]));
    #
    my $FA = $fields[6];
    $FA =~ s/\%//;
    $FA /= 100;
    my $VAQ = ".";

    my $FET = sprintf("%e", $fields[11]); # p value from fishers exact test

    ##placeholder/dummy, some will be corrected by downstream tool
    my $dbsnp_id = ".";
    my $qual = ".";
    my $filter = "PASS";
    my $info = ".";
    ##need SS check in here for somatic status to come out properly..
    my $format = "GT:GQ:DP:BQ:MQ:AD:FA:VAQ:FET";

    # If the variant called is N, just null out the GT and FET fields to minimize interference with cross-sample VCFs
    if ($var_allele_iub eq "N") {
        $GT = "./.";
        $FET = ".";
        $alt_alleles = ".";
        $AD = ".";
    }

    my $sample_string =join (":", ($GT, $GQ, $DP, $BQ, $MQ, $AD, $FA, $VAQ, $FET));
    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref_allele, $alt_alleles, $qual, $filter, $info, $format, $sample_string);
    return $vcf_line;
}

sub get_format_meta {
    my $self = shift;

    # Get all of the base FORMAT lines
    my @tags = $self->SUPER::get_format_meta; 

    my $fet = {MetaType => "FORMAT", ID => "FET", Number => 1, Type => "String", Description => "P-value from Fisher's Exact Test"};

    return (@tags, $fet);
}

