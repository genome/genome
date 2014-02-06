#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Model::SomaticVariation::Command::TestHelpers qw( create_test_objects run_test );

my $pkg = 'Genome::Model::Tools::Analysis::ResolveIgvReferenceName';
use_ok($pkg);

subtest "GRCh37-lite-build37" => sub {
    my $cmd = $pkg->execute(
        reference_name => 'GRCh37-lite-build37',
    );
    ok($cmd, 'Command ran without errors');
    is($cmd->igv_reference_name, 'b37', 'IGV reference name correct');
};

subtest "test reference name" => sub {
    my $cmd = $pkg->execute(
        reference_name => 'test reference name',
    );
    ok($cmd, 'Command ran without errors');
    is($cmd->igv_reference_name, 'unknown', 'IGV reference name correct');
};

done_testing();
