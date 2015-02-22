#!/usr/bin/env perl

use above 'Genome';

use Test::More;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Build;
use Genome::Utility::Test qw/compare_ok/;

use strict;
use warnings;

my $pkg = 'Genome::Model::RnaSeq::Command::ErccSpikeInQc';

use_ok($pkg);

# path = /gscmnt/gc13003/info/test_suite_data/Genome-Model-RnaSeq-Command-ErccSpikeInQc/v1
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "v1");

my $ercc_spike_in_file = $data_dir .'/ERCC_Controls_Analysis.txt';
my $expected_output_dir = $data_dir .'/expected_output'; 
my $build_data_dir = $data_dir .'/build';

my $rnaseq_model = Genome::Test::Factory::Model::RnaSeq->setup_object();
my $rnaseq_build = Genome::Test::Factory::Build->setup_object(
    model_id        => $rnaseq_model->id,
    data_directory  => $build_data_dir,
    status          => "Succeeded",
);

my $output_dir = Genome::Sys->create_temp_directory();  
my $cmd = $pkg->create(
   models => $rnaseq_model,
   ercc_spike_in_mix => 1,
   ercc_spike_in_file => $ercc_spike_in_file,
   output_directory => $output_dir,
);
ok($cmd,'create '. $pkg );
ok($cmd->execute(),'execute '. $pkg );       

compare_ok($expected_output_dir .'/summary.tsv', $output_dir .'/'. $rnaseq_model->id .'/summary.tsv');   

done_testing();