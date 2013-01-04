#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use YAML;
use Test::More;
use File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 5;
}

use_ok('Genome::Model::Tools::Sam::InferPairedStatus');

my $data_dir = __FILE__ . '.d';

my $paired_cmd = Genome::Model::Tools::Sam::InferPairedStatus->create(input=>$data_dir . "/paired.bam");
ok($paired_cmd->execute(), "executed cmd on data we expect to be paired");
ok($paired_cmd->is_paired_end, "is successfully paired like we expected it");

my $fragment_cmd = Genome::Model::Tools::Sam::InferPairedStatus->create(input=>$data_dir . "/fragment.bam");
ok($fragment_cmd->execute(), "executed cmd on data we expect to be fragment");
ok(! $fragment_cmd->is_paired_end, "is successfully fragment like we expected it");


done_testing();
