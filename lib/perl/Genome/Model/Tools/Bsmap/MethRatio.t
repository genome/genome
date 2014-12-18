#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;

my $class = 'Genome::Model::Tools::Bsmap::MethRatio';
use_ok($class);

sub create_test_object {
    my $temp_output_dir = Genome::Sys->create_temp_directory();
    my $example_bam = File::Spec->join($temp_output_dir, 'foo.bam');
    return $class->create(
        bam_file => $example_bam,
        output_directory => $temp_output_dir,
        reference => '36',
    );
}

{
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
}
