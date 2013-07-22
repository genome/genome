#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;
use Genome::Command::Tester qw(run_and_diff);

run_and_diff(
    command => '$script_dir/converge/convergeDrugGenes.pl '
        . ' --model_group_id=31779 --event_types_list=all  '
        . ' --dgidb_subdir_name=drugbank  --filter_name=antineo  --outdir=$output_dir' 
        . ' --reference_annotations_dir=$annotation_dir  --gene_groups=LUC17 --verbose=1',
    results_version => '2013-02-28',
    eventual_class => 'Genome::Model::ClinSeq::Command::Converge::DrugGenes',
);

