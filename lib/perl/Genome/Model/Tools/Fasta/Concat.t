#!/usr/bin/env genome-perl

use warnings;
use strict;

use above 'Genome';
use Test::More tests => 6;
use File::Temp qw/tempdir/;

my $class = 'Genome::Model::Tools::Fasta::Concat';
use_ok($class);

my %fasta_data = (
    file1 => {
        seq1_1 => "A" x80,
        seq1_2 => "AC"x40,
    },
    file2 => {
        seq2_1 => "C" x80,
        seq2_2 => "CG"x40,
    },
    file3 => {
        seq3_1 => "G" x80,
        seq3_2 => "GT"x40,
    },
    file4 => {
        seq4_1 => "T" x80,
        seq4_2 => "TA"x40,
    },
);

my $tmpdir = tempdir(CLEANUP => 1);
my $expected_path = "$tmpdir/expected.fa";
my $actual_path = "$tmpdir/actual.fa";

# create input files
my @input_files;
my $expected_fh = new IO::File(">$expected_path") || 
    die "Failed to open $expected_path for writing.";

for my $file (sort keys %fasta_data) {
    my $path = "$tmpdir/$file.fa";
    my $fh = new IO::File(">$path") || die "Failed to open $path for writing.";
    for my $seq (sort keys %{$fasta_data{$file}}) {
        $fh->print(">$seq\n$fasta_data{$file}{$seq}\n\n");
        $expected_fh->print(">$seq\n$fasta_data{$file}{$seq}\n");
    }
    $fh->close();
    push(@input_files, $path);
}
$expected_fh->close();

my $cmd = $class->create(input_files => \@input_files, output_file => $actual_path);
ok($cmd, 'Created command');
ok($cmd->execute, 'Executed command');
my $diff_result = `diff $expected_path $actual_path 2>&1`;
is($diff_result, '', 'Result is as expected');

$cmd = $class->create(input_files => [$input_files[0], $input_files[0]], output_file => $actual_path);
ok($cmd, 'Created command');
eval { $cmd->execute; };
ok($@, 'Trying to merge files with duplicate sequence ids is an error');

done_testing();
