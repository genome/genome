#!/usr/bin/env perl

use above 'Genome';

use Data::Dumper;
use Genome::File::Vcf::Entry;
use Genome::File::Vcf::Reader;
use File::Basename qw/basename/;
use Genome::Utility::Vcf 'diff_vcf_file_vs_file';
use Test::More;
use File::Slurp qw/read_file/;

use strict;
use warnings;

my $pkg = "Genome::Model::Tools::Joinx::VcfAnnotateMulti";
use_ok($pkg);

my $test_data_directory = __FILE__ . '.d';

my $input_file = join("/", $test_data_directory, "in.clean.vcf");
my $dbsnp_vcf = join("/", $test_data_directory, "dbsnp.clean.vcf");
my $thousand_genomes_vcf = join("/", $test_data_directory, "1kg-wgs.clean.vcf");
my $expected_file = join("/", $test_data_directory, "out.clean.vcf");

my $temp_dir = Genome::Sys->create_temp_directory();
my $output_file = join("/", $temp_dir, "output.vcf.gz");

my @specs = (
    Genome::Model::Tools::Joinx::VcfAnnotationSpec->create(
        annotation_file => $dbsnp_vcf,
        identifiers => 1,
        info_fields => ["dbSNPBuildID=dbsnp,per-alt", "GMAF"],
        ),
    Genome::Model::Tools::Joinx::VcfAnnotationSpec->create(
        annotation_file => $thousand_genomes_vcf,
        identifiers => 0,
        info_fields => ["AA", "AF", "AMR_AF", "ASN_AF", "EUR_AF", "VT"],
        ),
);

my $bad_cmd = $pkg->create(annotation_specs => \@specs, input_file => $input_file,
    bgzip_output => 1);

eval {
    $bad_cmd->execute;
};
ok($@, "Specifying bgzip_output without output file is an error");

my $cmd = $pkg->create(annotation_specs => \@specs, input_file => $input_file,
    output_file => $output_file, bgzip_output => 1);
$cmd->execute();
ok(-s $output_file, "Output file exists and is not empty");
ok(Genome::Sys->file_is_gzipped($output_file), "Output file is compressed");

my $reader = Genome::File::Vcf::Reader->new($output_file);
ok($reader, "Created vcf reader");
my $header = $reader->header;
ok($header, "Got vcf header");


# Strip all but the filename from the annotation line header
# The current version of joinx adds a trailing comma (probably a typo),
# but I plan to take it out in the next release.
my @annotation_files =
    sort
    map { $_ =~ s/^##annotation=.*\/([^,]*),?/$1/; $_ }
    grep { /^##annotation/ }
    $header->lines;

my @expected_annotation_files = sort map {basename($_->annotation_file)} @specs;
is_deeply(\@annotation_files, \@expected_annotation_files, "Annotation sources are listed in vcf header");

my $diff = diff_vcf_file_vs_file($output_file, $expected_file, ignore_patterns => ["^##annotation="]);
ok(!$diff, "Output matches expected value") or diag("Diff results:\n$diff");

done_testing();
