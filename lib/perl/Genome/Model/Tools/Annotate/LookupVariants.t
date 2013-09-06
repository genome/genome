#!/usr/bin/env genome-perl

use strict;
use warnings;
use File::Temp;
use above "Genome";
use Test::More tests => 12;

use_ok('Genome::Model::Tools::Annotate::LookupVariants');
my $variant_file = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-LookupVariants/snp.list.in";
ok (-e $variant_file);

my $tmpdir = File::Temp::tempdir(
    TEMPLATE => 'Model-Tools-Annotate-LookupVariants-XXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
);

my $dbsnp_path = "/gscmnt/sata835/info/medseq/model_data/2857166586/ImportedVariations/tmp";
my $known_out = "$tmpdir/Genome-Model-Tools-Annotate-LookupVariants-known-only.out";
if ($known_out && -e $known_out) {`rm -f $known_out`;}
my $exp_known_out = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-LookupVariants/expected.known-only.out";
ok (-e $exp_known_out);

my $known = Genome::Model::Tools::Annotate::LookupVariants->create(report_mode => "known-only", variant_file => "$variant_file", output_file => "$tmpdir/Genome-Model-Tools-Annotate-LookupVariants-known-only.out", filter_out_submitters => "SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG", require_allele_match => 1, dbSNP_path => $dbsnp_path);
ok ($known);

my $kv = $known->execute;

ok ($kv);
ok (-e $known_out);
my $knowndiff = `diff $exp_known_out $known_out`;
ok($knowndiff eq '', "known output as expected");


my $novel_out ="$tmpdir/Genome-Model-Tools-Annotate-LookupVariants-novel-only.out";
if ($novel_out && -e $novel_out) {`rm -f $novel_out`;}
my $exp_novel_out = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-LookupVariants/expected.novel-only.out";
ok (-e $exp_novel_out);


my $novel = Genome::Model::Tools::Annotate::LookupVariants->create(report_mode => "novel-only", variant_file => "$variant_file", output_file => "$tmpdir/Genome-Model-Tools-Annotate-LookupVariants-novel-only.out", filter_out_submitters => "SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG", require_allele_match => 1, dbSNP_path => $dbsnp_path);
ok ($novel);

my $nv = $novel->execute;
ok ($nv);
ok (-e $novel_out);
my $noveldiff = `diff $exp_novel_out $novel_out`;
ok($noveldiff eq '', "novel output as expected");
