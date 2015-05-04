#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Sys;
use File::Spec;
use File::Temp;
use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::UpdateFileForImport') or die;
use_ok('Genome::Model::Tools::CgHub::Test') or die;
Genome::Model::Tools::CgHub::Test->overload_lwp_user_agent_request;

my $legacy_sample_id = 'TCGA-L5-A8NH-01A-11R-A37I-31';
my $analysis_id = '01f22763-6bb2-4edc-b65a-99904f7c6fad';
my %expected_data = (
    analysis_id => $analysis_id,
    'individual.name' => join('-', (split('-', $legacy_sample_id))[0..2]),
    'individual.taxon' => 'human',
    'individual.upn' => join('-', (split('-', $legacy_sample_id))[1..2]),
    'sample.name' => $legacy_sample_id,
    'sample.nomenclature' => 'TCGA',
    'sample.common_name' => 'tumor',
    'sample.extraction_type' => 'rna',
    'sample.tissue_desc' => 'esophagus', 
    'sample.disease_abbr' => 'ESCA', 
    'library.name' => $legacy_sample_id.'-extlibs',
    'instdata.analysis_id' => $analysis_id,
);
my @expected_headers = sort keys %expected_data;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);

my $file = File::Spec->join($tmpdir, 'input.csv');
my $fh = Genome::Sys->open_file_for_writing($file);
$fh->print("analysis_id,source_files\n");
$fh->print("$analysis_id,file.bam\n");
$fh->close;

my $output_file = File::Spec->join($tmpdir, 'output.csv');
my $cmd = Genome::Model::Tools::CgHub::UpdateFileForImport->execute(file => $file, output_file => $output_file);
ok($cmd->result, "execute with file => $file");
my $content = Genome::Sys->read_file($output_file);
is_deeply(
    $content, 
    join(",", @expected_headers, 'source_files')."\n".join(",", map({ $expected_data{$_} } @expected_headers), 'file.bam')."\n",
    "correct content for file",
);

# analysis_id not found
throws_ok(sub{ $cmd->_resolve_common_name('TCGA'); }, qr/No sample type for analysis id: TCGA/, '_resolve_common_name w/ unknown analysis_id');
throws_ok(sub{ $cmd->_resolve_extraction_type('TCGA'); }, qr/No analyte code for analysis id: TCGA/, '_resolve_extraction_type w/ unknown analysis_id');
is($cmd->_resolve_tissue_desc('TCGA'), 'unknown', '_resolve_tissue_desc w/ unknown analysis_id gives unknown');
$cmd->delete;

done_testing();
