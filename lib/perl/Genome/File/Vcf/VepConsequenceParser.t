#!/usr/bin/env perl

use above 'Genome';

use Data::Dumper;
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Header;
use Genome::File::Vcf::Entry;
use IO::String;
use Test::More;

use strict;
use warnings;

my $pkg = 'Genome::File::Vcf::VepConsequenceParser';
use_ok($pkg);

# This data is clean. It is taken from dbsnp with some made up stuff thrown in.
my $vcf_fh = new IO::String(<<EOS
##fileformat=VCFV4.0
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|HGNC|DISTANCE|SIFT|PolyPhen|CELL_TYPE|Condel">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
1	69588	rs141776804	GGT	TGA,GG,G,GC	.	.	CSQ=-|ENSG00000186092|ENST00000335137|Transcript|FRAMESHIFT_CODING&feature_truncation|499-500|499-500|167||||OR4F5|||||,GA|ENSG00000186092|ENST00000335137|Transcript|NON_SYNONYMOUS_CODING|499-500|499-500|167|V/D|GTc/GAc||OR4F5||deleterious(0)|possibly_damaging(0.871)||deleterious(0.783),C|ENSG00000186092|ENST00000335137|Transcript|FRAMESHIFT_CODING|499-500|499-500|167||||OR4F5|||||,G|ENSG00000186092|ENST00000335137|Transcript|FRAMESHIFT_CODING|499-500|499-500|167||||OR4F5|||||	.
1	10177	.	A	C	.	.	CSQ=C|ENSG00000223972|ENST00000456328|Transcript|UPSTREAM|||||||DDX11L1|1692||||,C|ENSG00000227232|ENST00000488147|Transcript|DOWNSTREAM|||||||WASH7P|4227||||,C|ENSG00000227232|ENST00000541675|Transcript|DOWNSTREAM|||||||WASH7P|4186||||,C|ENSG00000223972|ENST00000450305|Transcript|UPSTREAM|||||||DDX11L1|1833||||,C|ENSG00000223972|ENST00000515242|Transcript|UPSTREAM|||||||DDX11L1|1695||||
EOS
);

my $reader = Genome::File::Vcf::Reader->fhopen($vcf_fh, "Test Vcf");
ok($reader, "Created vcf reader");
my $header = $reader->header;
ok($header, "Got vcf header");
ok(exists $header->info_types->{CSQ}, "CSQ tag exists in header");

my $csq = new $pkg($header);
ok($csq, "Created consequence parser");
ok($csq->isa($pkg), "Consequence parser has correct class");

my $multi_alt = $reader->next;
ok($multi_alt, "Got multi-alt vcf entry");
my $multi_transcript = $reader->next;
ok($multi_transcript, "Got multi-transcript vcf entry");

subtest "multiple_alts" => sub {
    my $e = $multi_alt;
    my $rv = $csq->process_entry($e);
    ok($rv, "Parsed VEP consequence");
    my @alleles = sort keys %$rv;
    my @expected = sort @{$e->{alternate_alleles}};
    is_deeply(\@alleles, \@expected, "Saw expected alleles");

    # Check expected values for the 2bp deletion
    subtest "2bp deletion values" => sub {
        my $data = $rv->{G};
        is(1, @$data, "Found one record for 2bp deletion");
        $data = $data->[0];

        my %expected = (
            allele => "-",
            gene => "ENSG00000186092",
            feature => "ENST00000335137",
            feature_type => "Transcript",
            consequence => "FRAMESHIFT_CODING&feature_truncation",
            cdna_position => "499-500",
            cds_position => "499-500",
            protein_position => 167,
            amino_acids => '',
            codons => '',
            existing_variation => '',
            hgnc => "OR4F5",
            distance => '',
            sift => '',
            polyphen => '',
            cell_type => '',
            condel => '',
        );

        is_deeply($data, \%expected, "Expected values for 2bp deletion")
            or diag("Got: " . Dumper($data) . ", Expected: " . Dumper(\%expected));
    };

    subtest "1bp deletion values" => sub {
        my $data = $rv->{GG};
        is(1, @$data, "Found one record for 1bp deletion");
        $data = $data->[0];

        my %expected = (
            allele => "G",
            gene => "ENSG00000186092",
            feature => "ENST00000335137",
            feature_type => "Transcript",
            consequence => "FRAMESHIFT_CODING",
            cdna_position => "499-500",
            cds_position => "499-500",
            protein_position => 167,
            amino_acids => '',
            codons => '',
            existing_variation => '',
            hgnc => "OR4F5",
            distance => '',
            sift => '',
            polyphen => '',
            cell_type => '',
            condel => '',
        );

        is_deeply($data, \%expected, "Expected values for 2bp deletion")
            or diag("Got: " . Dumper($data) . ", Expected: " . Dumper(\%expected));
    };

    subtest "GGT -> TGA" => sub {
        my $data = $rv->{TGA};
        is(1, @$data, "Found one record for GGT -> TGC");
        $data = $data->[0];

        my %expected = (
            allele => "GA", # Verifying that VEP does the wrong thing.
            gene => "ENSG00000186092",
            feature => "ENST00000335137",
            feature_type => "Transcript",
            consequence => "NON_SYNONYMOUS_CODING",
            cdna_position => "499-500",
            cds_position => "499-500",
            protein_position => 167,
            amino_acids => 'V/D',
            codons => 'GTc/GAc',
            existing_variation => '',
            hgnc => "OR4F5",
            distance => '',
            sift => 'deleterious(0)',
            polyphen => 'possibly_damaging(0.871)',
            cell_type => '',
            condel => 'deleterious(0.783)',
        );

        is_deeply($data, \%expected, "Expected values for 2bp deletion")
            or diag("Got: " . Dumper($data) . ", Expected: " . Dumper(\%expected));
    };
};

subtest "multiple_transcripts" => sub {
    my $e = $multi_transcript;
    my $data = $csq->process_entry($e);
    my @alleles = keys %$data;
    is_deeply(\@alleles, ["C"], "allele matches");
    is(5, @{$data->{C}}, "Got 5 records");

    my @expected = (
        {
            allele => 'C',
            gene => 'ENSG00000223972', feature => 'ENST00000456328',
            feature_type => 'Transcript', consequence => 'UPSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "DDX11L1",
            distance => '1692',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
        {
            allele => 'C',
            gene => 'ENSG00000227232', feature => 'ENST00000488147',
            feature_type => 'Transcript', consequence => 'DOWNSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "WASH7P",
            distance => '4227',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
        {
            allele => 'C',
            gene => 'ENSG00000227232', feature => 'ENST00000541675',
            feature_type => 'Transcript', consequence => 'DOWNSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "WASH7P",
            distance => '4186',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
        {
            allele => 'C',
            gene => 'ENSG00000223972', feature => 'ENST00000450305',
            feature_type => 'Transcript', consequence => 'UPSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "DDX11L1",
            distance => '1833',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
        {
            allele => 'C',
            gene => 'ENSG00000223972', feature => 'ENST00000515242',
            feature_type => 'Transcript', consequence => 'UPSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "DDX11L1",
            distance => '1695',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
    );

    is_deeply($data->{C}, \@expected, "Got expected transcripts");
};

subtest "filters" => sub {
    my $distance_filter = sub {
        my $data = shift;
        return $data->{distance} > 4000;
    };

    my $filtering_parser = new $pkg($header, filters => [$distance_filter]);
    my $e = $multi_transcript;
    my $rv = $filtering_parser->process_entry($e);
    ok($rv, "Parsed VEP consequence");

    my @expected = (
        {
            allele => 'C',
            gene => 'ENSG00000227232', feature => 'ENST00000488147',
            feature_type => 'Transcript', consequence => 'DOWNSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "WASH7P",
            distance => '4227',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
        {
            allele => 'C',
            gene => 'ENSG00000227232', feature => 'ENST00000541675',
            feature_type => 'Transcript', consequence => 'DOWNSTREAM',
            cdna_position => '', cds_position => '', protein_position => '',
            amino_acids => '', codons => '',  existing_variation => '',
            hgnc => "WASH7P",
            distance => '4186',
            sift => '', polyphen => '', cell_type => '', condel => ''
        },
    );

    is_deeply($rv->{C}, \@expected, "Got expected transcripts");

};

subtest "format transcripts" => sub {
    my $e = $multi_transcript;
    my $rv = $csq->process_entry($e);

    my $fmt = $csq->format_transcripts($rv);
    my $raw_info = $e->info("CSQ");
    my @expected_fields = sort split(",", $raw_info);
    my @actual_fields = sort split(",", $fmt);
    is_deeply(\@actual_fields, \@expected_fields, "transcript reformatting");

    my $empty_fmt = $csq->format_transcripts({});
    is(".", $empty_fmt, "empty transcripts => .");
};

done_testing();
