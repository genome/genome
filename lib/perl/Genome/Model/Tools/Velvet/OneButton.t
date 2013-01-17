#!/usr/bin/env genome-perl
use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";
use Test::More;
use File::Basename qw(dirname basename);
use File::Compare qw(compare);

my $archos = `uname -a`;
if ($archos !~ /64/) {
    plan skip_all => "Must run from 64-bit machine";
}

# OLD:
# oneButtonVelvet-3opt.pl 1k-trimmed.fastq -i 260 -g 4500000 --hash 31 33 35 --version 0.7.57-64 -o old

# NEW:
# gmt velvet one-button 1k-trimmed.fastq -i 260 -g 4500000 --hash 31,33,35 --version 0.7.57-64 -o new

my $module = 'Genome::Model::Tools::Velvet::OneButton';
use_ok($module) or die;

# invalid version (checked first)
ok(!Genome::Model::Tools::Velvet::OneButton->create(version => 2), 'Failed to create as expectred w/ invalid version');

# CHANGE THIS WHENEVER WE INTENTIONALLY SWITCH OUTPUT
my $version = 'v3';

my @data_sub_dirs_orig = ("v1-a", "v1-b",);#"v1-c");

my @data_sub_dirs = ("$version-a", "$version-b",);#"$version-c");
my @params = ('-i 260 -g 4500000 --hash 31,33,35 --version 0.7.57-64 --min-contig-length 100 -o output-dir',
    '-i 260 -g 4500000 --hash 31,33,35  --bound-enumeration 2 --min-contig-length 100 --version 0.7.57-64 -o output-dir',);
  # '-i 260 -g 4500000 --hash 31,33,35 --c 17,19,21 --bound-enumeration 2 --version 0.7.57-64 -o output-dir');

#make sure @data_sub_dirs and @params contain same # of elements

foreach my $param (@params) {

    my $data_dir = $module;
    $data_dir =~ s/::/-/g;
    $data_dir = $ENV{GENOME_TEST_INPUTS} . "/$data_dir";

    #UNCOMMENT THE FOLLOWING TWO LINES WHEN BUILD NEW TEST DATA
    #`mkdir -p $data_dir/$data_sub_dirs[0]`;
    #`cp -rf $data_dir/$data_sub_dirs_orig[0]/input.fastq $data_dir/$data_sub_dirs[0]/.`;

    ok(-d $data_dir, "found data directory $data_dir");
    my $data_sub_dir = '/'.shift @data_sub_dirs;
    $data_dir .= $data_sub_dir;

    my $expected_dir = $data_dir . '/output-dir';
    ok(-d $expected_dir, "found expected data directory $expected_dir");

    my $expected_stdout = $data_dir . '/actual.stdout';
    ok(-e $expected_stdout, "found expected data directory $expected_stdout");

    my $expected_stderr = $data_dir . '/actual.stderr';
    ok(-e $expected_stderr, "found expected data directory $expected_stderr");

    my $input_file = "$data_dir/input.fastq";
    ok(-e $input_file, "found input file $input_file");

    my $temp_dir = Genome::Sys->create_temp_directory();
    # SWITCH TO THIS WHEN WE WANT TO GENERATE INTENTIONALLY NEW TEST DATA
    #my $temp_dir = "$ENV{PWD}/velvet$data_sub_dir";
    #`mkdir -p $temp_dir`;
    ok(-d $temp_dir, "temp directory made at $temp_dir");

    my $actual_dir = "$temp_dir/output-dir";
    mkdir $actual_dir;
    ok(-d $actual_dir, "made output dir $actual_dir");

    my $actual_stdout = "$temp_dir/actual.stdout";
    my $actual_stderr = "$temp_dir/actual.stderr";

    my $cmd = "chdir $temp_dir; gmt velvet one-button $input_file $param > actual.stdout 2>actual.stderr";
    note($cmd);
    my $rv = system($cmd);
    $rv /= 256;
    ok($rv == 0, "command runs successfully");

    for my $old_file (glob("$actual_dir/*input.fastq*")) {
    }

    for my $test_file (glob("$actual_dir/*")) {
        if ($test_file =~ /input\.fastq/) {
            $test_file = replace_date($test_file);
        }
        diff_test($test_file, $actual_dir, $expected_dir, $temp_dir);
    }

    my @stdout_diff = `sdiff -s $expected_stdout $actual_stdout | grep -v -- '$temp_dir' | grep -vP '^Reading FastQ file'`;
    is(scalar(@stdout_diff), 0, "stdout matches")
        or diag(@stdout_diff);

    my @stderr_diff = `sdiff -s $expected_stderr $actual_stderr | grep -v -- '$temp_dir' | grep -v GENOME_DEV_MODE`;
    is(scalar(@stderr_diff), 0, "stderr matches except for the line with a date and DEV_MODE notice")
        or diag(@stderr_diff);
}

done_testing();

sub replace_date {
    my $old_file = shift;
    my $dirname = dirname($old_file);
    my $basename = basename($old_file);
    $basename =~ s/^.*input.fastq/SOMEDATE-input.fastq/;
    my $new_file = "$dirname/$basename";
    rename $old_file, $new_file;
    ok(-e $new_file, "renamed $old_file to $new_file");
    return $new_file;
}

sub diff_test {
    my ($test_file, $actual_dir, $expected_dir, $temp_dir) = @_;

    (my $common_path = $test_file) =~ s/$actual_dir//;
    my $expected_file = join '/', $expected_dir, $common_path;

    if ($common_path =~ /Log$/) {
        my @velvetlog_diff = `sdiff -s $test_file $expected_file | grep -v -- '$temp_dir' | grep -v output-dir`;
        is(scalar(@velvetlog_diff), 2, "log diffed as expected: $common_path")
            or diag(@velvetlog_diff);
    } else {
        my $diff = sub {
            map { $_ =~ s/^#Bams: .*// } @_;
            map { $_ =~ s/^#input file: .*// } @_;
            map { $_ =~ s/.*velvet. output-dir .*// } @_;
            map { $_ =~ s/^\w{3} \w{3} \d\d \d\d:\d\d:\d\d \d{4}$// } @_;
            my $c = $_[0] ne $_[1];
            if ($c == 1) {
                diag("First diff:\n--- " . $_[0] . "+++ " . $_[1]);
            }
            return $c;
        };
        is(compare($test_file, $expected_file, $diff), 0, "output diffed as expected: $common_path") or diag qx(diff -u $test_file $expected_file);
    }
}
