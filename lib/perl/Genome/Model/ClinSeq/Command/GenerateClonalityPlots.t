#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
#use Test::More skip_all => "very slow";
use Test::More tests => 3; 

my $expected_out = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-ClinSeq-Command-GenerateClonalityPlots/2012-11-21.full';

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

ok(-d $expected_out, "directory of expected output exists: $expected_out");

# Run the tool as described in the synopsis.
my $cmd = "genome model clin-seq generate-clonality-plots --somatic-var-build=129396826  --output-dir=$actual_out  --common-name='AML54'  --verbose";
#$cmd .= " --limit=1000";
#$cmd .= " --read-counts=$expected_out/allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts"; 
eval { Genome::Sys->shellcmd(cmd => $cmd); };
ok(!$@, "executed command correctly: $cmd")
    or do {
        die("exception: $@");
    };

# The last column in the tsv file has output which varies randomly from run-to-run. :(
# We replace that value with ? before doing a diff so that we won't get spurious failures.
# In shell: cat AML54.clustered.data.tsv | perl -nae '$F[-1] = "?"; print join("\t",@F),"\n"' >| AML54.clustered.data.tsv.testmasked
my $fhin = Genome::Sys->open_file_for_reading("$actual_out/AML54.clustered.data.tsv");
my $fhout = Genome::Sys->open_file_for_writing("$actual_out/AML54.clustered.data.tsv.testmasked");
while (my $row = <$fhin>) {
    chomp $row;
    my @fields = split("\t",$row);
    $fields[-1] = "?";
    $fhout->print(join("\t",@fields),"\n");
}
$fhout->close;

# The differences test excludes files which always differ (embed dates, or are the subject of a masking as above).
my @diff = `diff -X $expected_out/IGNORE_IN_DIFFS -r --brief $expected_out $actual_out`;
is(scalar(@diff), 0, "no differences") 
    or do {
        for (@diff) { diag($_) }
        warn "*** if the above differences are not in error, rebuild the test data by running this test with REBUILD on the command-line ***";
    }

