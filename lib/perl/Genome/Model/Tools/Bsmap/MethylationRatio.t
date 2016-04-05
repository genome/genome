#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Slurp;
use File::Compare;
use List::AllUtils qw(max);
use Genome::Utility::Test;
use Test::More tests => 3;

my $class = 'Genome::Model::Tools::Bsmap::MethylationRatio';
use_ok($class);

my $test_root = Genome::Utility::Test->data_dir($class, 'v1');

sub fasta_for_reference {
    my ($name, $version) = @_;
    return Genome::Model::ImportedReferenceSequence->get(
        name => $name
    )->build_by_version($version)->full_consensus_path('fa');
}

sub create_test_object {
    my $temp_output_dir = Genome::Sys->create_temp_directory();
    
    my $bam_file = File::Spec->join($test_root, 'input', 'all_sequences.bam');
    my $reference = fasta_for_reference('TEST-human', '1');

    return $class->create(
        bam_file => $bam_file,
        output_directory => $temp_output_dir,
        reference => $reference
    );
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
    plan tests => 6;

    my $test_obj = create_test_object();

    is($test_obj->output_file, 'methylation', 'has correct default output_file');


    my $expected_path = File::Spec->join(
        $test_root, 'expected', 'methylation.bgz');
     my $expected_path_index = File::Spec->join(
        $test_root, 'expected', 'methylation.bgz.tbi');
    
    my $actual_path = File::Spec->join(
        $test_obj->output_directory, 'methylation.bgz');
    my $actual_path_index = File::Spec->join(
        $test_obj->output_directory, 'methylation.bgz.tbi');
        

    ok($test_obj->execute(), 'Command executes');

    ok(-s $actual_path, "Output $actual_path exists");
    
    ok(-s $actual_path, "Output $actual_path_index exists");
    is(compare($expected_path, $actual_path), 0,
        'Expected and derived files are equal');
    
    is(compare($expected_path_index, $actual_path_index), 0,
        'Expected and derived files indices are equal');
};
