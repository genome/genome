#!/usr/bin/env perl

use above 'Genome';
use Test::More;
use Genome::File::Vcf::Reader;
use Genome::File::BamReadcount::Entry;

use strict;
use warnings FATAL => 'all';


my $pkg = "Genome::File::Vcf::BamReadcountUtilities";
use_ok($pkg);

subtest "readcount position conversion" => sub {
    my $entry1 = create_vcf_entry(20, 1234567, 'GTC', 'G,ATC,GTCT'); 
    my $expected_positions = [1234567, 1234568, 1234569];
    is_deeply( [Genome::File::Vcf::BamReadcountUtilities::vcf_entry_to_readcount_positions($entry1)],
        $expected_positions,
        "Positions match with multiple alleles as expected");

    my $entry2 = create_vcf_entry(20, 1234567, 'G', 'A'); 
    my ($pos2) = Genome::File::Vcf::BamReadcountUtilities::vcf_entry_to_readcount_positions($entry2);
    is($pos2, 1234567, "SNP position matches");
    
    my $entry3 = create_vcf_entry(20, 1234567, 'GTC', 'G'); 
    my ($pos3) = Genome::File::Vcf::BamReadcountUtilities::vcf_entry_to_readcount_positions($entry3);
    is($pos3, 1234568, "Deletion position matches");
    
    my $entry4 = create_vcf_entry(20, 1234567, 'GTC', 'GTCT'); 
    my ($pos4) = Genome::File::Vcf::BamReadcountUtilities::vcf_entry_to_readcount_positions($entry4);
    is($pos4, 1234569, "Insertion position matches");
};

subtest "entry comparison" => sub {
    my $vcf_entry = create_vcf_entry(20, 1234567, 'GTC', 'G,ATC,GTCT'); 
    
    my $brct_entry1 = create_bam_readcount_entry(21, 1234569);
    ok(!Genome::File::Vcf::BamReadcountUtilities::entries_match($brct_entry1, $vcf_entry), "Don't match if on different chromosomes");
    
    my $brct_entry2 = create_bam_readcount_entry(20, 1234569);
    ok(Genome::File::Vcf::BamReadcountUtilities::entries_match($brct_entry2, $vcf_entry), "Match if one allele has same position");
    
    my $brct_entry3 = create_bam_readcount_entry(20, 1234550);
    ok(!Genome::File::Vcf::BamReadcountUtilities::entries_match($brct_entry3, $vcf_entry), "Don't match if no alleles have the same position");
};

sub create_vcf_entry {
    my ($chrom, $pos, $ref, $variant) = @_;

    my $vcf_str = <<EOS;
##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=noident,Description="No identifier">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
$chrom\t$pos\tmicrosat1,foo\t$ref\t$variant\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
EOS
    my $vcf_fh = new IO::String($vcf_str);
    my $reader = Genome::File::Vcf::Reader->fhopen($vcf_fh, "Test Vcf");

    my $vcf_entry = $reader->next;
    return $vcf_entry;
}

sub create_bam_readcount_entry {
    my ($chrom, $pos) = @_;

    my $brct_str = <<BRCT;
$chrom	$pos	T	98	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:1:42.00:2.00:0.00:0:1:0.78:0.04:14.00:1:0.13:201.00:0.13	C:1:30.00:2.00:0.00:1:0:0.55:0.13:42.00:0:nan:199.00:0.22	G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	T:96:36.73:26.54:0.00:50:46:0.49:0.02:17.53:92:0.41:229.41:0.40	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
BRCT
    my $brct_entry = Genome::File::BamReadcount::Entry->new($brct_str);
    return $brct_entry;
}


done_testing();
