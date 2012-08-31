#!/usr/bin/env genome-perl

use above "Genome";
use warnings;
use strict;
use Data::Dumper;

my @solr_docs;

my @models = Genome::Model->get();

my @models2 = Genome::Model->get(name => {operator=>'LIKE', value=>'%incremental%'});

my %output_to_cache;

for (@models) {

    my %store_hash;

    my $last_build = $_->last_complete_build;
    my $last_build_data_dir;
    my $last_build_id;
    if ($last_build) {

         $last_build_data_dir = $last_build->data_directory;
         $last_build_id = $last_build->id;
    }


    my $model_data_dir = $_->data_directory;

    my @builds = $_->builds;

    my $metric_count;
    for (@builds) {  
        $metric_count++ if ($_->get_metric('gold-heterozygous-snp match heterozygous-1-allele-variant filtered'));
    }   

    $store_hash{'gold_snp_het_metric_count'} = $metric_count;
    $store_hash{'last_complete_build_id'} = $last_build_id;
    $store_hash{'last_complete_build_data_directory'} = $last_build_data_dir;
    $store_hash{'data_directory'} = $model_data_dir;
    $store_hash{'subject_name'} = $_->subject_name;
    $store_hash{'user_name'} = $_->user_name;
    $store_hash{'creation_date'} = $_->creation_date;

    $output_to_cache{'search_model_metrics_' . $_->id} = \%store_hash;
}

print Dumper(\%output_to_cache);
