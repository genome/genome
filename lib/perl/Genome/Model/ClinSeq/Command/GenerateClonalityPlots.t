#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
#use Test::More skip_all => "very slow"; 
use Test::More tests => 4; 

my $expected_out = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-ClinSeq-Command-GenerateClonalityPlots/2013-04-04/';

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
my $somvar_build = Genome::Model::Build->get($somvar_build_id);
ok($somvar_build, "Got somatic variation build from id: $somvar_build_id") or die;

my $cmd = Genome::Model::ClinSeq::Command::GenerateClonalityPlots->create(somatic_var_build=>$somvar_build, output_dir=>$actual_out, common_name=>'HCC1395', verbose=>1, chromosome=>'22');
$cmd->queue_status_messages(1);
my $r1 = $cmd->execute();
is($r1, 1, 'Testing for successful execution.  Expecting 1.  Got: '.$r1);


# The last column in the tsv file has output which varies randomly from run-to-run. :(
# We replace that value with ? before doing a diff so that we won't get spurious failures.
# In shell: cat AML54.clustered.data.tsv | perl -nae '$F[-1] = "?"; print join("\t",@F),"\n"' >| AML54.clustered.data.tsv.testmasked
#my $fhin = Genome::Sys->open_file_for_reading("$actual_out/AML54.clustered.data.tsv");
#my $fhout = Genome::Sys->open_file_for_writing("$actual_out/AML54.clustered.data.tsv.testmasked");
#while (my $row = <$fhin>) {
#    chomp $row;
#    my @fields = split("\t",$row);
#    $fields[-1] = "?";
#    $fhout->print(join("\t",@fields),"\n");
#}
#$fhout->close;

# The differences test excludes files which always differ (embed dates, or are the subject of a masking as above).
my @diff = `diff -x '*.pdf' -r $expected_out $actual_out`;
is(scalar(@diff), 11, "only expected differences") 
or do {
  for (@diff) { diag($_) }
  warn "*** if the above differences are not in error, rebuild the test data by running this test with REBUILD on the command-line ***";
  my $temp_dir = "/tmp/last-generate-clonality-plots-result/";
  Genome::Sys->shellcmd(cmd => "rm -fr $temp_dir");
  Genome::Sys->shellcmd(cmd => "mv $actual_out $temp_dir");
}


