#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;

use above "Genome";
use Cache::Memcached;

my $lock_resource = $ENV{GENOME_LOCK_DIR} . '/gcsearch/memcache_loader';

my $lock = Genome::Sys->lock_resource(resource_lock=>$lock_resource, max_try=>0);
unless ($lock) {
    die "could not lock, another instance must be running.";
}


# 3 hours, this'll run every two via cron
my $expiration = 60*60*3;

main();
Genome::Sys->unlock_resource(resource_lock=>$lock);
exit;


sub main {

    my $memd = new Cache::Memcached {
        'servers' => ['lims-dev:11211'],
        'debug'   => 0,
    };

    my $genome_stuff = get_genome_stuff();
    for my $key (keys %$genome_stuff) {
        $memd->set($key, $genome_stuff->{$key}, $expiration);
    }

    my $cn = get_common_names();
    $memd->set('genome_individual', $cn, $expiration);

}

sub get_common_names {

    return {
        AML3  => [ 2816877266, 2816874825 ],
        AML4  => [ 2816901765, 2816901764 ],
        AML9  => [ 2816902483, 2816905138 ],
        AML10 => [ 2816905223, 2816903316 ],
        AML11 => [ 2816909095, 2816908646 ],
        AML12 => [ 2816909096, 2816912119 ],
    };
}

sub get_genome_stuff {

    my @solr_docs;

    my @pps = Genome::ProcessingProfile->get();
    my @models = Genome::Model->get();

    my @models2 = Genome::Model->get(name => {operator=>'LIKE', value=>'%incremental%'});

    my %output_to_cache;

    for (@models) {

        my %store_hash;

        eval {
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
        };
    }

    return \%output_to_cache;

}





