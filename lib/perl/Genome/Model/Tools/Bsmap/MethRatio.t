#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Slurp;
use List::AllUtils qw(max);
use Test::More tests => 3;

my $class = 'Genome::Model::Tools::Bsmap::MethRatio';
use_ok($class);

my $test_root = File::Spec->join($ENV{GENOME_TEST_INPUTS},
    qw(Genome-Model-Tools-Bsmap-MethRatio v2.74 2015-01-06));

sub fasta_for_reference {
    my ($name, $version) = @_;
    return Genome::Model::ImportedReferenceSequence->get(
        name => $name
    )->build_by_version($version)->fasta_file();
}

sub create_test_object {
    my $temp_output_dir = Genome::Sys->create_temp_directory();
    my $bam_file = File::Spec->join($test_root, 'input', 'all_sequences.bam');
    my $reference = fasta_for_reference('TEST-human', '1');

    return $class->create(
        bam_file => $bam_file,
        output_directory => $temp_output_dir,
        reference => $reference,
        version => '2.74',
    );
}

sub diff_paths {
    my ($expected_path, $actual_path) = @_;

    my @expected = read_file($expected_path);
    my @actual = read_file($actual_path);

    my $lines = max($#expected, $#actual);
    my $diffs = 0;
    for (0 .. $lines) {
        if ($actual[$_] ne $expected[$_]) {
            $diffs++;
        }
    }

    return $diffs;
}

subtest 'Dereference reference sequence' => sub {
    plan tests => 2;

    my $test_obj = create_test_object();

    my %mapping = (
        36 => 'NCBI-human-build36',
        37 => 'GRCh37-lite-build37'
    );

    for my $ref_short_name (keys %mapping) {
        $test_obj->reference($ref_short_name);

        my $fasta = Genome::Model::Build::ReferenceSequence->get(
            name => $mapping{$ref_short_name}
        )->cached_full_consensus_path('fa');

        is($test_obj->_reference_fasta, $fasta, 'Got correct fasta object');
    }
};

subtest 'Execute' => sub {
    plan tests => 4;

    my $test_obj = create_test_object();
    my $test_output_name = 'test.snvs.hq';

    is($test_obj->output_file, 'snvs.hq', 'has correct default output_file');

    $test_obj->output_file($test_output_name);

    my $expected_path = File::Spec->join(
        $test_root, 'expected', $test_output_name);

    my $actual_path = File::Spec->join(
        $test_obj->output_directory, $test_output_name);

    ok($test_obj->execute(), 'Command executes');

    ok(-s $actual_path, "Output $actual_path exists");

    is(diff_paths($expected_path, $actual_path), 0,
        'Found no diffs between expected and actual test.snvs.hq');
};
