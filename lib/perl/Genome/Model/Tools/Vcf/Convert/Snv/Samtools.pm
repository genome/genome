package Genome::Model::Tools::Vcf::Convert::Snv::Samtools;

use strict;
use warnings;
use Genome;
use POSIX 'strftime';
use Genome::Info::IUB;

class Genome::Model::Tools::Vcf::Convert::Snv::Samtools {
    is  => 'Genome::Model::Tools::Vcf::Convert::Base',
    doc => 'Generate a VCF file from samtools output',
};

sub help_synopsis {
    <<'HELP';
    Generate a VCF file from samtools snv output
HELP
}

sub help_detail {
    <<'HELP';
    Parses the input file and creates a VCF containing all the snvs.
HELP
}

sub source {
    return 'Samtools';
}

#single sample for now
sub _get_header_columns {
    my $self = shift;
    my @header_columns = ("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",$self->aligned_reads_sample);
    return @header_columns;
}


sub parse_line {
    my ($self, $line) = @_;
    return if $line =~ /^#/; # no mpileup vcf header here
    my @columns = split("\t", $line);

    if ($columns[4] =~ /^[A-Z]+/) { #mpileup output vcf format already
        my @alt = split /,/, $columns[4];
        pop @alt if $alt[$#alt] eq 'X';  #remove the silly X in ALT colum if it exists

        #my ($dp, $fa, $dp4, $mq) = $columns[7] =~ /DP=(\d+);AF1=(\S+?);\S+DP4=(\S+?);MQ=(\d+);/;
        my @info = split /;/, $columns[7];
        my %info;
        for my $info (@info) {
            my ($key, $value) = split /=/, $info;
            $info{$key} = $value;
        }

        my ($dp, $mq, $dp4, $af1) = ($info{DP}, $info{MQ}, $info{DP4}, $info{AF1});
        unless (defined $dp and defined $mq and defined $dp4 and defined $af1) {
            die $self->error_message("Fail to get valid info from line: $line");
        }

        my @sample = split /:/, $columns[9];
        my ($gt, $gq) = ($sample[0], $sample[-1]);

        $columns[4] = join ',', @alt;
        $columns[6] = 'PASS';
        $columns[7] = '.';
        $columns[8] = 'GT:GQ:DP:MQ:DP4:FA:BQ:SS'; #no way to calculate BQ, VAQ(snp quality)
        $columns[9] = join ':', $gt, $gq, $dp, $mq, $dp4, $af1, '.', '.';

        my $col_ct   = scalar $self->_get_header_columns;
        my $new_line = join "\t", splice(@columns, 0, $col_ct); #remove some unwanted columns in test
        return $new_line;
    }

    my ($chr, $pos, $ref, $genotype, $gq, $vaq, $mq, $dp, $read_bases, $base_quality) = @columns;
    #replace ambiguous/IUPAC bases with N in ref
    $ref =~ s/[^ACGTN\-]/N/g;

    my @alt_alleles = Genome::Info::IUB->variant_alleles_for_iub($ref, $genotype);
    my @alleles     = Genome::Info::IUB->iub_to_alleles($genotype);
    my @all_alleles = (uc($ref), @alt_alleles);

    my $alt = join ",", @alt_alleles;

    #add the ref and alt alleles' positions in the allele array to the GT field
    my $gt = $self->generate_gt($ref, \@alt_alleles, \@alleles);

    # Parse the pileup and quality strings so that they have the same length and can be mapped to one another
    if ($read_bases =~ m/[\$\^\+-]/) {
        $read_bases =~ s/\^.//g; #removing the start of the read segement mark
        $read_bases =~ s/\$//g; #removing end of the read segment mark
        while (my ($match, $indel_len) = $read_bases =~ m/([\+-]{1}(\d+))/g) {
            my $pos = index($read_bases, $match);
            substr($read_bases, $pos, $indel_len + length($match), "");
            #$read_bases =~ s/[\+-]{1}$indel_len.{$indel_len}//; # remove indel info from read base field
        }
    }

    my (%ad, %bq_total);
    my ($bq_string, $ad_string);
    if ( length($read_bases) != length($base_quality) ) {
        die $self->error_message("After processing, read base string and base quality string do not have identical lengths: $read_bases $base_quality");
    }

    # Count the number of times the variant occurs in the pileup string (AD) and its quality from the quality string (BQ)
    my @bases     = split("", $read_bases);
    my @qualities = split("", $base_quality);
    my $total_ad  = 0;

    for (my $index = 0; $index < scalar(@bases); $index++) {
        my $base = uc($bases[$index]);
        $base    = uc($ref) if $base eq '.' or $base eq ',';
        for my $variant (@all_alleles) {
            if ($variant eq $base) {
                $ad{$variant}++;
                $total_ad++ unless $variant eq uc($ref);
                #http://samtools.sourceforge.net/pileup.shtml base quality is the same as mapping quality
                $bq_total{$variant} += ord($qualities[$index]) - 33; 
            }
        }
    }

    # Get an average of the quality for BQ
    my %bq;
    
    for my $variant (@all_alleles) {
        if ($ad{$variant}) {
            $bq{$variant} = int($bq_total{$variant} / $ad{$variant});
        } 
        else {
            $bq{$variant} = 0;
            $ad{$variant} = 0;
        }
    }
    $bq_string = join ",", map { $bq{$_} } @all_alleles;
    $ad_string = join ",", map { $ad{$_} } @all_alleles;

    # fraction of reads supporting alt
    my $fa = $total_ad / $dp; 
    $fa = sprintf "%.3f", $fa; # Round to 3 decimal places since we dont have that many significant digits

    # If the variant called is N, just null out the GT and FET fields to minimize interference with cross-sample VCFs
    if ($genotype eq "N") {
        $gt  = "./.";
        $alt = "N";
        $ad_string = ".";
        $bq_string = ".";
        $fa = ".";
    }

    # Placeholder for later adjustment
    my $dbsnp_id = ".";
    my $qual = $vaq;
    my $filter = "PASS";
    my $format = "GT:GQ:DP:BQ:MQ:AD:FA:VAQ:SS";
    my $info = ".";
    my $ss   = ".";
    my $sample_string = join (":", ($gt, $gq, $dp, $bq_string, $mq, $ad_string, $fa, $vaq, $ss));

    my $vcf_line = join("\t", $chr, $pos, $dbsnp_id, $ref, $alt, $qual, $filter, $info, $format, $sample_string);

    return $vcf_line;
}

