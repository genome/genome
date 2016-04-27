#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use above 'Workflow';

use Test::More tests => 5;
use File::Compare;
use File::Copy;

BEGIN {
        $ENV{NO_LSF} = 1;
        use_ok ('Genome::Model::Tools::Blat::Subjects');
}
my $data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Blat-Subjects';
my $query_file = $data_dir .'/test.fa';
my $expected_psl = $data_dir .'/test.psl';
my $tmp_dir = Genome::Sys->create_temp_directory;

my $tmp_query_file = $tmp_dir . '/test.fa';
copy $query_file, $tmp_query_file;

my $psl_path = $tmp_dir .'/test_tmp.psl';
my $blat_output_path = $tmp_dir .'/test_tmp.out';
my $blat_params = '-mask=lower -out=pslx -noHead';
my $ref_seq_dir = Genome::Config::reference_sequence_directory() . '/refseq-for-test';

opendir(DIR,$ref_seq_dir) || die "Failed to open dir $ref_seq_dir";
my @ref_seq_files = map { $ref_seq_dir .'/'. $_ } grep { !/^all_seq/ } grep { /\.fa$/ } readdir(DIR);
closedir(DIR);

is(scalar(@ref_seq_files),3,'expected three input subject files');

my $blat = Genome::Model::Tools::Blat::Subjects->create(
                                                        query_file => $tmp_query_file,
                                                        subject_files => \@ref_seq_files,
                                                        psl_path => $psl_path,
                                                        blat_params => $blat_params,
                                                        blat_output_path => $blat_output_path
                                                  );
isa_ok($blat,'Genome::Model::Tools::Blat::Subjects');
ok($blat->execute,'execute '. $blat->command_name);
ok(!compare($psl_path,$expected_psl),'psl files are the same');
