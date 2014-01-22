#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use Genome::Command::Tester qw(run_and_diff);

run_and_diff(
    command => '$script_dir/converge/convergeSnvs.pl  --model_group_id=66909 --outdir=$output_dir  --label=BRAF  --cancer_annotation_db=tgi/cancer-annotation/human/build37-20131010.1  --verbose=1',
    results_version => '2013-03-06',
    eventual_class => 'Genome::Model::ClinSeq::Command::Converge::Snvs',
);

