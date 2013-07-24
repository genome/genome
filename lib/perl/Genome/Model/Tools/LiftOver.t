#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Genome::Command::Tester qw(run_and_diff);
use Test::More tests => 13;

use_ok('Genome::Model::Tools::LiftOver');

run_and_diff(
    command => 'gmt lift-over --input-is-annoformat --source $input_dir/snvs.anno --dest $output_dir/snvs.hg18ToHg19.liftover.anno --unmapped $output_dir/snvs.hg18ToHg19.lost.anno --lift-direction hg18ToHg19',
    results_version => '2013-07-19-snv',
);

run_and_diff(
    command => 'gmt lift-over --input-is-annoformat --source $input_dir/indels.anno --dest $output_dir/indels.hg18ToHg19.liftover.anno --unmapped $output_dir/indels.hg18ToHg19.lost.anno --lift-direction hg18ToHg19',
    results_version => '2013-07-19-indel',
);

run_and_diff(
    command => 'gmt lift-over --file-format svold --source $input_dir/svs.svold --dest $output_dir/svs.hg18ToHg19.liftover.svold -unmapped $output_dir/svs.hg18ToHg19.lost.svold --lift-direction hg18ToHg19',
    results_version => '2013-07-19-sv-old',
);

