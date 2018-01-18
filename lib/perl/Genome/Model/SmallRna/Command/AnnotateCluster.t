#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Utility::Test;
use File::Spec;

my $class = 'Genome::Model::SmallRna::Command::AnnotateCluster';
use_ok($class);

my $data_dir    = Genome::Utility::Test->data_dir_ok($class);
my $test_file   = Genome::Sys->create_temp_file_path;
my $exp_out_dir = $data_dir . '/v1';

my $cluster_bed_file = $data_dir.'/top_sorted_clusters.bed';
my $exp_out_file     = $exp_out_dir.'/annotation_intersect.tsv';

my $annot_files = join(',', map {
        File::Spec->join( $data_dir, 'Annotation_Files', $_)
    }
    qw(Curated_clusters_FINAL_ANNOTATION.bed miRBasev16_BioMart_annotation_combined.bed NCBI-human.combined-annotation_58_37c_v2_modified.bed)
);

my $annot_name = 'Curated_Clusters,miRBase_BioMart,TGI';

my $cmd = $class->execute(
    cluster_bed_file    => $cluster_bed_file,
    annotation_bed_file => $annot_files,
    annotation_name     => $annot_name,
    output_tsv_file     => $test_file,
);

Genome::Utility::Test::compare_ok($test_file, $exp_out_file, 'output file annotation_intersect.tsv created ok');

done_testing;

