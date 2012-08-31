#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 5;
use Test::MockObject;
use File::Temp qw/ tempdir /;
use File::Compare;
use Data::Dumper;

# Check test data directory exists
my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ImportTranscriptFeature-KeggPathway';
ok (-d $test_dir, "Test data directory exists");

# Check test gene pathway file exists
my $test_gene_file = $test_dir . '/gene_pathway_test.tsv';
ok (-e $test_gene_file, "Test gene pathway file exists");

# Check test path description file exists
my $test_path_file = $test_dir . '/pathway_desc_test.tsv';
ok (-e $test_path_file, "Test pathway description file exists");

# Mock the kegg wsdl
my ($server) = Test::MockObject->new;
$server->fake_module('SOAP::Lite', service => sub {return $server});

my @mock_paths = ( {entry_id => 0, definition => 'mock path'} );
$server->set_always('list_pathways', \@mock_paths);

my @mock_genes = qw/ hsa:00 hsa:01 /;
$server->set_always('get_genes_by_pathway', \@mock_genes);

# Create and execute kegg pathway import tool
my $output_dir = "$ENV{GENOME_TEST_TEMP}";
my $temp_dir = tempdir('/Genome-Model-Tools-ImportTranscriptFeature-KeggPathway-XXXXXX',
    DIR => $output_dir, CLEANUP => 1);
my $gene_output = $temp_dir . '/gene_path.tsv';
my $path_output = $temp_dir . '/path_desc.tsv';

my $kegg_obj = Genome::Model::Tools::ImportTranscriptFeature::KeggPathway->create(
    gene_pathway_file => $gene_output,
    pathway_description_file => $path_output,
);

$kegg_obj->execute;

ok (compare($test_gene_file, $gene_output) == 0, "Gene/pathway output matches test output");
ok (compare($test_path_file, $path_output) == 0, "Pathway/pathway description output matches test output");

unlink $gene_output;
unlink $path_output;
