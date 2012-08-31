package Genome::Model::Tools::Vcf::Convert::Indel::Varscan;

use strict;
use warnings;
use Genome;
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Indel::Varscan {
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
    my $leading_base = $fields[2];
    my $dbsnp_id = ".";
    my $indel_string = $fields[18]; #last field on varscan line should be: +ACT or -GA deletion/insertion reporting style
    my $alt_alleles;
    my $ref_allele;
    if($indel_string =~ m/\+/)  { #insertion
     $alt_alleles = $leading_base . substr($indel_string, 1);
     $ref_allele = $leading_base;
    }
    elsif($indel_string =~m/-/) { #deletion
        #we should switch ref and alt (compared to insertions) to show that those bases are going away
        $ref_allele = $leading_base . substr($indel_string,1);
        $alt_alleles = $leading_base;

    }
    else {
        $self->error_message("line does not *appear* to contain a +/- char in final field, don't know what to do: $line");
        return 0;
    }
    #TODO this is turned off for now because it interferes with applying filters (bed coordinates will be different once left shifted)
    # ($chr, $pos, $ref_allele, $alt_alleles) = $self->normalize_indel_location($chr, $pos, $ref_allele, $alt_alleles);
    my $genotype_call = $fields[3];
    my $GT;
    if($genotype_call =~ m/\*/) { #my belief is that varscan only calls het or hom. so if we see match an asterisk, we know its het. 
        $GT="0/1";
    }else {
        $GT="1/1";
    }
    my $DP=$fields[4]+$fields[5];  #ref allele supporting reads + var allele supporting reads
    my $MQ= '.'; #$fields[13]; # ref allele supporitng read map qual
    my $AD = $fields[5];
    my $FA = sprintf("%0.3f",$fields[5]/$DP);
    my $FET = sprintf("%e", $fields[11]);


    ##fill in defaults for some shit
    my $BQ=".";
    my $GQ=".";
    my $VAQ = ".";
    my $filter = "PASS";
    my $info = ".";
    my $qual = ".";

    ##need SS check in here for somatic status to come out properly..
    my $format = "GT:GQ:DP:BQ:MQ:AD:FA:VAQ:FET";
    my $sample_string =join (":", ($GT, $GQ, $DP, $BQ, $MQ, $AD, $FA, $VAQ, $FET));
    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref_allele, $alt_alleles, $qual, $filter, $info, $format, $sample_string);
    return $vcf_line;
}

sub get_format_meta {
    my $self = shift;

    # Get all of the base FORMAT lines
    my @tags = ($self->common_format_meta, $self->extra_format_meta); 

    my $fet = {MetaType => "FORMAT", ID => "FET", Number => 1, Type => "String", Description => "P-value from Fisher's Exact Test"};

    return (@tags, $fet);
}




