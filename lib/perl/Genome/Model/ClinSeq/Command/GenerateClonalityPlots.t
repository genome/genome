#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 8;

my $expected_out =
    Genome::Config::get('test_inputs') . '/Genome-Model-ClinSeq-Command-GenerateClonalityPlots/2014-08-28/';

# "REBUILD" on the command-line when running the test sets the output to be the reference data set
# IMPORTANT: change the directory above to a new date when updating test data!
my $actual_out;
if (@ARGV and $ARGV[0] eq 'REBUILD') {
    warn "**** rebuilding expected output at $expected_out ****";
    mkdir $expected_out unless -d $expected_out;
    $actual_out = $expected_out;
}
else {
    $actual_out = Genome::Sys->create_temp_directory;
}

ok(-d $expected_out, "directory of expected output exists: $expected_out") or die;

#Get a somatic variation build
my $somvar_build_id = 135798051;
my $somvar_build    = Genome::Model::Build->get($somvar_build_id);
ok($somvar_build, "Got somatic variation build from id: $somvar_build_id") or die;

my $cmd = Genome::Model::ClinSeq::Command::GenerateClonalityPlots->create(
    bam_readcount_version => 0.6,
    somatic_var_build     => $somvar_build,
    misc_annotation_db    => Genome::Db->get("tgi/misc-annotation/human/build37-20130113.1"),
    chromosome            => '22',
    verbose               => 1,
    output_dir            => $actual_out,
    common_name           => 'HCC1395',
);
$cmd->queue_status_messages(1);
my $r1 = $cmd->execute();
is($r1, 1, 'Testing for successful execution.');

#Since we can not diff the pdf files, at least check for file creation...
for my $pdf_name (
    qw(
    HCC1395.clonality.pdf
    HCC1395.clonality.cn2.pdf
    HCC1395.clonality.filtered_snvs.pdf
    HCC1395.clonality.filtered_snvs.cn2.pdf
    )
    )
{
    my $pdf_path = join('/', $actual_out, $pdf_name);
    ok(-s $pdf_path, "Found non-zero PDF file $pdf_name");
}

# The differences test excludes files which always differ (embed dates, or are the subject of a masking as above).
my $temp_dir = "/tmp/last-generate-clonality-plots-result/";
my @diff     = `diff -x '*.pdf' -x '*.R' -r $expected_out $actual_out`;
is(scalar(@diff), 0, "only expected differences")
    or do {
    for (@diff) {diag($_)}
    warn
        "*** if the above differences are not in error, rebuild the test data by running this test with REBUILD on the command-line after updating the expected dir's date ***";
    Genome::Sys->shellcmd(cmd => "rm -fr $temp_dir");
    Genome::Sys->shellcmd(cmd => "mv $actual_out $temp_dir");
    };

