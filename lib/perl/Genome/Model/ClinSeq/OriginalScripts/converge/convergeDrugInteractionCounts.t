#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use Genome::Command::Tester qw(run_and_diff);

#TODO: switch to 65642 ?

run_and_diff(
    command => '$script_dir/converge/convergeDrugInteractionCounts.pl '
        . ' --model_group_id=44083 --event_types_list=all  --dgidb_subdir_names=drugbank,santa_monica_lung '
        . ' --filter_name=default  --outdir=$output_dir  --verbose=1',
    results_version => '2013-02-28',
    eventual_class => 'Genome::Model::ClinSeq::Command::Converge::DrugInteractionCounts',
);

