#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw();
use File::Basename qw(dirname);
use File::Slurp qw(read_file);
use File::chdir qw($CWD);

my $pkg = "Genome::Model::Tools::Bwa::RunMem";
use_ok($pkg);

my $TEST_BWA_VERSION = "0.7.10";

## Test data setup
my $temp_dir = Genome::Sys->create_temp_directory;
my $data_dir = Genome::Utility::Test->data_dir($pkg, "v1");
my @input_fastqs = glob("$data_dir/*.fq");
is(scalar @input_fastqs, 2, "Found 2 input files");
my $header_file = sprintf("%s/header.sam", $data_dir);
ok(-s $header_file, "Sam header file exists");
my $fasta_file = sprintf("%s/ref/ref.fa", $data_dir);
ok(-s $fasta_file, "Test reference sequence exists");

my @rg_lines = grep {/^\@RG/} read_file($header_file);
is(scalar @rg_lines, 1, "Found read group in header");
my $read_group_line = $rg_lines[0];
chomp $read_group_line;
$read_group_line =~ s/\t/\\t/g;

my $samtools = $pkg->_samtools_path;

## Object params that do not vary between tests
my %common_params = (
    aligner_params => "-R '$read_group_line'",
    indexed_fasta => $fasta_file,
    aligner_index_fasta => $fasta_file,
    sam_header_path => $header_file,
    max_sort_memory => "200M",
    num_threads => 1,
    samtools_version => '0.1.19',
    );

subtest "invalid samtools versions" => sub {
    my $obj = $pkg->create(
            bwa_version => "0.5.9",
            input_fastqs => \@input_fastqs,
            output_file => "/dev/null",
            aligner_log_path => "/dev/null",
            %common_params
        );

    ok($obj, "created command");

    my @invalid = ("r123", "0.1.18", "1.0");

    for my $ver (@invalid) {
        $obj->samtools_version($ver);
        eval { $obj->_validate_params; };
        ok($@, "samtools_version = [$ver] is an error");
    }

    $obj->samtools_version("0.1.19");
    eval { $obj->_validate_params; };
    ok(!$@, "samtools version 0.1.19 is ok");
};

subtest "specifying threads with -t is an error" => sub {
    my $obj = $pkg->create(
            bwa_version => "0.5.9",
            input_fastqs => \@input_fastqs,
            output_file => "/dev/null",
            aligner_log_path => "/dev/null",
            %common_params
        );

    ok($obj, "created command");

    my @invalid = ("-t 4", "-t4", "-S -t4", "-S -t 4 -e",
        "-t '4'", '-t "4"', '-t " 4"', "-t ' 4  '");
    my @valid = ("-R '\@RG\tID:rg-test'");
    for my $params (@invalid) {
        $obj->aligner_params($params);
        eval { $obj->_validate_params; };
        ok($@, "aligner params == [$params] is an error");
    }

    for my $params (@valid) {
        $obj->aligner_params($params);
        eval { $obj->_validate_params; };
        ok(!$@, "aligner params == [$params] is not an error");
    }
};

subtest "bwa versions without mem" => sub {
    my $obj = $pkg->create(
            bwa_version => "0.5.9",
            input_fastqs => \@input_fastqs,
            output_file => "/dev/null",
            aligner_log_path => "/dev/null",
            %common_params
        );

    ok($obj, "created command");
    eval { $obj->execute; };
    ok($@, "specifying bwa versions that lack mem is an error");
};

subtest "too many input files" => sub {
    my $obj = $pkg->create(
            bwa_version => $TEST_BWA_VERSION,
            input_fastqs => [qw(1 2 3)],
            output_file => "/dev/null",
            aligner_log_path => "/dev/null",
            %common_params
        );

    ok($obj, "created command");
    eval { $obj->execute; };
    ok($@, "too many input files is an error");
};


subtest "paired-end alignment" => sub {
    my $expected_pe = sprintf("%s/expected-pe.bam", $data_dir);
    my $output_file = sprintf("%s/pe.bam", $temp_dir);
    my $log_file = sprintf("%s/pe.log", $temp_dir);

    my $obj = $pkg->create(
            bwa_version => $TEST_BWA_VERSION,
            input_fastqs => \@input_fastqs,
            output_file => $output_file,
            aligner_log_path => $log_file,
            %common_params
        );

    ok($obj, "created command");
    ok($obj->execute, "Executed command");
    ok(-s $output_file, "Output file is not empty");
    ok(-s $log_file, "Log file is not empty");

    compare_sam($output_file, $expected_pe);
};

subtest "single-end alignment" => sub {
    my $expected_se = sprintf("%s/expected-se.bam", $data_dir);
    my $output_file = sprintf("%s/se.bam", $temp_dir);
    my $log_file = sprintf("%s/se.log", $temp_dir);

    my $obj = $pkg->create(
            bwa_version => $TEST_BWA_VERSION,
            input_fastqs => [grep {/r1\.fq/} @input_fastqs],
            output_file => $output_file,
            aligner_log_path => $log_file,
            %common_params
        );

    ok($obj, "created command");
    ok($obj->execute, "Executed command");
    ok(-s $output_file, "Output file is not empty");
    ok(-s $log_file, "Log file is not empty");

    compare_sam($output_file, $expected_se);
};

## Failure tests
#
# There was some concern expressed that we might miss failures when
# some part of the pipeline fails in this command. We're explicitly
# setting pipefail in bash, so this shouldn't be an issue. Let's
# test to be sure!
#
# We'll patch each method that generates a command in the pipeline
# with a call to perl -e 'exit 1;' and test that we get an exception.

my $perl_interp = $^X;
my $fail_cmd = "$perl_interp -e 'exit 1;'";

sub make_failure_test {
    my ($name, $method) = @_;
    return sub {
        # We don't want to leave core files laying around
        local $CWD = $temp_dir;

        no strict 'refs';
        local *$method = sub { return $fail_cmd; };
        use strict;

        my $output_file = sprintf("%s/failure.bam", $temp_dir);
        my $log_file = sprintf("%s/failure.log", $temp_dir);

        my $obj = $pkg->create(
                bwa_version => $TEST_BWA_VERSION,
                input_fastqs => \@input_fastqs,
                output_file => $output_file,
                aligner_log_path => $log_file,
                %common_params
            );

        ok($obj->$method =~ /exit/, "bad things are going to happen");
        ok($obj, "created command");

        eval { $obj->execute; };
        ok($@, "Command failed due to $name crashing");
    };
}

my %failures = (
    "bwa" => "${pkg}::_aligner_command",
    "header replacement" => "${pkg}::_sam_replace_header_cmdline",
    "bam conversion" => "${pkg}::_sam_to_uncompressed_bam_cmdline",
    "sort" => "${pkg}::_sort_cmdline",
    "callmd" => "${pkg}::_calmd_cmdline",
    );

while (my ($name, $method) = each %failures) {
    subtest "$name failure detection" => make_failure_test($name, $method);
}

# Helper to compare raw sam output for equivalence
sub compare_sam {
    # A is the output file, B is expected
    my ($path_a, $path_b) = @_;

    my @lines_a = qx{$samtools view -h $path_a};
    my @lines_b = qx{$samtools view -h $path_b};
    chomp @lines_a;
    chomp @lines_b;

    is(scalar @lines_a, scalar @lines_b, "Same number of lines in sam files");

    my $last_chr;
    my $last_pos = -1;
    my $ordered = 1;
    my $ln = 0;
    for my $a (@lines_a) {
        ++$ln;
        next if $a =~ /^\@/;

        my @fields = split("\t", $a);
        my ($chr, $pos) = @fields[(2, 3)];
        if (!defined $last_chr or $last_chr ne $chr) {
            $last_chr = $chr;
            $last_pos = -1;
        }
        if ($pos < $last_pos) {
            $ordered = 0;
            diag "In $path_a:$ln: Bad sort order ".
                "($last_chr:$last_pos comes before $chr:$pos";
            last;
        }
    }
    ok($ordered, "sort order");

    # Now check that we have the same data.
    # If the SAM sorting program changes (e.g., from stable to unstable sort),
    # the order could be different (but equivalent) so we sort here.
    @lines_a = sort @lines_a;
    @lines_b = sort @lines_b;
    is_deeply(\@lines_a, \@lines_b, "Sam files contain the same data");

    return 1;
}

done_testing();
