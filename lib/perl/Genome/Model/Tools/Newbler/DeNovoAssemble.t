#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

my $archos = `uname -a`;
unless ($archos =~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

use_ok ( 'Genome::Model::Tools::Newbler::DeNovoAssemble' ) or die;

my $version = 'v2';
my $test_suite_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/DeNovoAssemble-'.$version;
ok ( -d $test_suite_dir, 'Test suite dir exists' ) or die;

my $input_fastq = $test_suite_dir.'/2869511846-input.fastq';

ok ( -s $input_fastq, 'Input fastq file exists' ) or die;

my $temp_dir = Genome::Sys->create_temp_directory();
ok ( -d $temp_dir, "Created temp test dir at $temp_dir" ) or die;

my $create = Genome::Model::Tools::Newbler::DeNovoAssemble->create(
    version          => 'mapasm454_source_03152011',
    output_directory => $temp_dir,
    input_files      => [ $input_fastq ],
    consed           => 1,
    rip              => 1,
    );
ok ( $create, "Created tool" ) or die;

ok ( $create->execute, "Executed tool" ) or die;

#files that sould match exactly
my @files_to_compare = qw/
    454AlignmentInfo.tsv  454AllContigs.qual  454LargeContigs.fna	454ReadStatus.txt
    454AllContigs.fna     454ContigGraph.txt  454LargeContigs.qual     	454TrimStatus.txt
/;
for my $file ( @files_to_compare ) {
    ok ( File::Compare::compare( $temp_dir."/$file", $test_suite_dir."/$file" ) == 0, "$file files match" );
}

#files that don't match due to time stamp and temp dir locations
my %files_with_diffs = (
    '454NewblerMetrics.txt' => 4,
    '454NewblerProgress.txt' => 2,
);
for my $file (sort keys %files_with_diffs) {
    my @diffs = `sdiff -s $temp_dir/$file $test_suite_dir/$file`;
    my $expected_num_diffs = $files_with_diffs{$file};
    ok ( scalar @diffs == $expected_num_diffs, "Correctly found 2 differences in $file" );
}

#ace and phdball files should match bec default time stamp used
for my $file ( '/consed/edit_dir/454Contigs.ace.1', '/consed/phdball_dir/phd.ball.1' ) {
    ok ( File::Compare::compare( $temp_dir."$file", $test_suite_dir."$file" ) == 0, "$file files match" );
}

#<STDIN>;

done_testing();
