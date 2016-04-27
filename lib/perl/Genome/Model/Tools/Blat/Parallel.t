#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use above 'Workflow';

use Test::More tests => 5;
#use Test::More skip_all => 'workflow and lsf issues taking a long time to test this';
use File::Compare;
use File::Temp;

BEGIN {
    use_ok ('Genome::Model::Tools::Blat::Parallel');
}

Genome::Config::set_env('workflow_builder_backend', 'inline');

my $tmp_dir = File::Temp::tempdir('Genome-Model-Tools-Blat-Parallel-XXXXX', CLEANUP => 1, TMPDIR => 1);
my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Blat-Parallel';

my @query_files = (
    $data_dir .'/s_7_1_sequence.fa',
    $data_dir .'/s_7_2_sequence.fa',
);
my $expected_psl = $data_dir .'/test.psl';
my $psl_path = $tmp_dir .'/test_tmp.psl';
my $blat_output_path = $tmp_dir .'/test_tmp.out';

my $blat_params = '-mask=lower -out=pslx -noHead';
my $ref_seq_dir = Genome::Config::reference_sequence_directory() . '/refseq-for-test';

opendir(DIR,$ref_seq_dir) || die "Failed to open dir $ref_seq_dir";
my @ref_seq_files = map { $ref_seq_dir .'/'. $_ } grep { !/^all_seq/ } grep { /\.fa$/ } readdir(DIR);
closedir(DIR);

is(scalar(@ref_seq_files),3,'expected three input subject files');

my $blat = Genome::Model::Tools::Blat::Parallel->create(
    query_files => \@query_files,
    subject_files => \@ref_seq_files,
    output_directory => $tmp_dir,
    psl_path => $psl_path,
    blat_params => $blat_params,
    blat_output_path => $blat_output_path,
);
isa_ok($blat,'Genome::Model::Tools::Blat::Parallel');
ok($blat->execute,'execute '. $blat->command_name);
ok(!compare($psl_path,$expected_psl),'psl files are the same');
