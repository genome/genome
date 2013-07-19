#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Genome::Command::Tester qw(run_and_diff);
use Test::More tests => 9;

use_ok('Genome::Model::Tools::LiftOver');

run_and_diff(
    command => 'gmt lift-over --input-is-annoformat --source $input_dir/snvs.anno --dest $output_dir/snvs.hg18ToHg19.liftover.anno --unmapped $output_dir/snvs.hg18ToHg19.lost.anno --lift-direction hg18ToHg19',
    results_version => '2013-07-19-snv',
);

run_and_diff(
    command => 'gmt lift-over --input-is-annoformat --source $input_dir/indels.anno --dest $output_dir/indels.hg18ToHg19.liftover.anno --unmapped $output_dir/indels.hg18ToHg19.lost.anno --lift-direction hg18ToHg19',
    results_version => '2013-07-19-indel',
);

#$ARGV[0] = ('REBUILD');
#run_and_diff(
#    command => 'gmt lift-over --file-format svpair --source $input_dir/svs.svpair --dest $output_dir/svs.hg18ToHg19.liftover.svpair -unmapped $output_dir/svs.hg18ToHg19.lost.svpair --lift-direction hg18ToHg19',
#    results_version => '2013-07-19-sv',
#);

