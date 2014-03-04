#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::GenotypeMicroarray::GenotypeFile::EntryFactory') or die;
use_ok('Genome::File::Vcf::Entry') or die;
use_ok('Genome::File::Vcf::Header') or die;

my ($vcf_entry);

my $entry_factory = Genome::Model::GenotypeMicroarray::GenotypeFile::EntryFactory->create(
    sample_name => '__TEST_SAMPLE__',
);
my $entry_from_hash = $entry_factory->build_entry(genotype_hash());
is_deeply($entry_from_hash, expected_entry_from_hash(), 'got entry from hash');

my $entry_from_vcf = $entry_factory->build_entry(vcf_entry());
is_deeply($entry_from_vcf, expected_entry_from_vcf(), 'got entry from vcf');

done_testing();

###

sub genotype_hash {
    return {
        'alleles' => 'AG',
        'log_r_ratio' => '-0.3639',
        'position' => '752566',
        'cnv_confidence' => 'NA',
        'allele2' => 'G',
        'reference' => 'A',
        'cnv_value' => '2.0',
        'chromosome' => '1',
        'sample_id' => '2879594813',
        'allele1' => 'A',
        'id' => 'rs3094315',
        'gc_score' => '0.8931'
    }
}

sub expected_entry_from_hash {
    return bless( 
        {
            'reference_allele' => 'A',
            'log_r_ratio' => '-0.3639',
            'position' => '752566',
            'chrom' => '1',
            'identifiers' => [ 'rs3094315' ],
            'cnv_confidence' => 'NA',
            'reference' => 'A',
            'cnv_value' => '2.0',
            'chromosome' => '1',
            '_format' => [ 'GT', 'OG' ],
            'sample_name' => '__TEST_SAMPLE__',
            'allele1' => 'A',
            'info_fields' => {
                'hash' => {
                    'LR' => '-0.3639',
                    'CV' => '2.0',
                    'GC' => '0.8931'
                },
                'order' => [ 'CV', 'LR', 'GC' ]
            },
            'id' => 'rs3094315',
            'quality' => '.',
            'alternate_alleles' => [ 'G' ],
            'gc_score' => '0.8931',
            '_sample_data' => [ [ '0/1', 'AG' ] ],
            'alleles' => 'AG',
            'allele2' => 'G',
            '_filter' => [],
            'sample_id' => '2879594813',
            '_format_key_to_idx' => {
                'OG' => 1,
                'GT' => 0
            }
        }, 'Genome::File::Vcf::Entry' );
}

sub vcf_entry {
    return $vcf_entry if $vcf_entry;
    my @lines = (
        '##fileformat=VCFv4.1',
        '##INFO=<ID=CC,Number=1,Type=Float,Description="CNV Confidence">',
        '##INFO=<ID=LR,Number=1,Type=Float,Description="Log R Ratio">',
        '##INFO=<ID=CV,Number=1,Type=Float,Description="CNV Value">',
        '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency for each ALT allele in the same order as listed">',
        '##INFO=<ID=GC,Number=1,Type=Float,Description="GC Score">',
        '##FORMAT=<ID=AL,Number=1,Type=String,Description="Alleles">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	__TEST_SAMPLE__'
    );
    my $header = Genome::File::Vcf::Header->create(lines => \@lines);
    my $line =  join("\t", (qw| 1 752566 rs3094315 A G . . CV=2.0;LR=-0.3639;GC=0.8931 GT:OG 0/1:AG |));
    $vcf_entry = Genome::File::Vcf::Entry->new($header, $line);
    return $vcf_entry;
}

sub expected_entry_from_vcf {
    my $entry = vcf_entry();
    $entry->info;
    $entry->sample_data;
    $entry->{sample_name} = '__TEST_SAMPLE__';
    $entry->{id} = $entry->{identifiers}->[0];
    $entry->{chromosome} = $entry->{chrom};
    $entry->{reference} = $entry->{reference_allele};
    my %info = (
        'allele_frequency' => 'NA',
        'allele1' => 'A',
        'allele2' => 'G',
        'alleles' => 'AG',
        'cnv_confidence' => 'NA',
        'cnv_value' => '2.0',
        'gc_score' => '0.8931',
        'log_r_ratio' => '-0.3639',
    );
    for my $k ( keys %info ) {
        $entry->{$k} = $info{$k};
    }

    return $entry;
}

