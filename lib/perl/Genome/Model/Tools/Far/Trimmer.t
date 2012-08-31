#!/usr/bin/env genome-perl

use strict;
use warnings;
use above "Genome";  # forces a 'use lib' when run directly from the cmdline
use Test::More tests => 15;
use File::Temp;

my $datadir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Far-Trimmer';
 
my $source1 = "$datadir/1.fastq";
my $source2 = "$datadir/2.fastq";
die "no sources" unless -e $source1 and -e $source2;

my $temp_dir = Genome::Sys->create_temp_directory;

#first test single adapter removal
my $cmd = Genome::Model::Tools::Far::Trimmer->create(
    source => $source1,
    source2 => $source2,
    target => "$temp_dir/single_adapter_removed.out",
    adapter => 'gtttcccagtcacgata',
    trim_end => 'any',
    cut_off => 2,
    min_overlap => 16,
    max_uncalled => 100,
    nr_threads => 4,
    format => 'fastq',
    use_version => '2.0',
    far_output => '/dev/null',
);
my @expected_output = ("$datadir/single_adapter_removed_1.fastq","$datadir/single_adapter_removed_2.fastq");

ok($cmd, "successfully created dust fastq command");
ok($cmd->execute, "successfully executed dust fastq command");
my @output = glob("$temp_dir/*fastq");
is(@output,2,"got 2 fastq files");
is(Genome::Sys->md5sum($output[0]), Genome::Sys->md5sum($expected_output[0]), 'first fastq file matches output');
is(Genome::Sys->md5sum($output[1]), Genome::Sys->md5sum($expected_output[1]), 'first fastq file matches output');

#now test adapter and reverse complement removal
my $cmd2 = Genome::Model::Tools::Far::Trimmer->create(
    source => $source1,
    source2 => $source2,
    target => "$temp_dir/single_adapter_removed.out",
    adapter => 'gtttcccagtcacgata',
    trim_reverse_complement => 1, 
    trim_end => 'any',
    cut_off => 2,
    min_overlap => 16,
    max_uncalled => 100,
    nr_threads => 4,
    use_version => '2.0',
    far_output => '/dev/null',
);
my @expected_output2 = ("$datadir/revcomp_adapter_removed_1.fastq","$datadir/revcomp_adapter_removed_2.fastq");

ok($cmd2, "successfully created dust fastq command");
ok($cmd2->execute, "successfully executed dust fastq command");
my @output2 = glob("$temp_dir/*fastq");
is(@output2,2,"got 2 fastq files");
is(Genome::Sys->md5sum($output2[0]), Genome::Sys->md5sum($expected_output2[0]), 'first fastq file matches output');
is(Genome::Sys->md5sum($output2[1]), Genome::Sys->md5sum($expected_output2[1]), 'first fastq file matches output');

#now test adapter and reverse complement removal with all params specified in params property
my $cmd3 = Genome::Model::Tools::Far::Trimmer->create(
    source => $source1,
    source2 => $source2,
    target => "$temp_dir/single_adapter_removed.out",
    params => '--trim-end any --cut-off 2 --adapter GTTTCCCAGTCACGATA --min-overlap 16 --max-uncalled 100 --nr-threads 4 --format fastq',
    trim_reverse_complement => 1, 
    use_version => '2.0',
    far_output => '/dev/null',
);
my @expected_output3 = ("$datadir/revcomp_adapter_removed_1.fastq","$datadir/revcomp_adapter_removed_2.fastq");

ok($cmd3, "successfully created dust fastq command");
ok($cmd3->execute, "successfully executed dust fastq command");
my @output3 = glob("$temp_dir/*fastq");
is(@output3,2,"got 2 fastq files");
is(Genome::Sys->md5sum($output3[0]), Genome::Sys->md5sum($expected_output3[0]), 'first fastq file matches output');
is(Genome::Sys->md5sum($output3[1]), Genome::Sys->md5sum($expected_output3[1]), 'first fastq file matches output');
