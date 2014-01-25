#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use File::Slurp qw(write_file);
use File::Temp qw(tempdir);
use JSON;
use Test::More;

use strict;
use warnings;

my $pkg = 'Genome::Model::PhenotypeCorrelation::Command::BeginBurdenTest';

use_ok($pkg);

my $tmpdir = tempdir(CLEANUP => 1);

my $glm_model_data = <<EOS
analysis_type	clinical_data_trait_name	variant/gene_name	covariates	memo
Q	quant1	NA	CV1+CV2+CV3	
B	cat1	NA	CV4+CV5+CV6	
EOS
;
my $glm_model_file_path = join("/", $tmpdir, "glm-model.txt");
write_file($glm_model_file_path, $glm_model_data);

my @genes = map {"GENE$_"} 1..2;
my %params = (
        genes => \@genes,
        option_file => "dummy.txt",
        glm_model_file => $glm_model_file_path,
        );

my $expected = [
    {
        gene => "GENE1",
        analysis_data_type => "Q",
        phenotype => "quant1",
        covariates => [map {"CV$_"} 1..3],
    },
    {
        gene => "GENE2",
        analysis_data_type => "Q",
        phenotype => "quant1",
        covariates => [map {"CV$_"} 1..3],
    },
    {
        gene => "GENE1",
        analysis_data_type => "B",
        phenotype => "cat1",
        covariates => [map {"CV$_"} 4..6],
    },
    {
        gene => "GENE2",
        analysis_data_type => "B",
        phenotype => "cat1",
        covariates => [map {"CV$_"} 4..6],
    },

];

my $cmd = $pkg->create(%params);
ok($cmd, "Created command");
ok($cmd->execute, "Executed command");

my $json = new JSON;
my @params = map {$json->decode($_)} $cmd->job_params;
is_deeply(\@params, $expected, "Params are as expected") or diag(
    "Expected: " . Dumper($expected) . "\n" .
    "Actual: " . Dumper(\@params));

done_testing();


