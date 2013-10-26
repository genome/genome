#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;
use File::Temp;
use File::Compare;

BEGIN {
    if (`uname -a` =~ /x86_64/){
        plan tests => 13;
    }
    else{
        plan skip_all => 'Must run on a 64 bit machine';
    }
    use_ok('Genome::Model::Tools::Sam::VarFilter');
}

my $root_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam/VarFilter';

my $tmp_dir  = File::Temp::tempdir(
    "VarFilter_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my $pileup_var_file = "$root_dir/test.pileup";
my $vcf_var_file    = "$root_dir/test.vcf";

my @out_files = qw(
    snv.varfilter
    indel.varfilter
    snv.varfilter.filtered
    indel.varfilter.filtered
);

my @formats = qw(pileup vcf);

my @pileup_out_files = map{$tmp_dir.'/'.$_.'.pileup'}@out_files;
my @vcf_out_files    = map{$tmp_dir.'/'.$_.'.vcf'}@out_files;

#using the default for samtools r963 samtools.pl varFilter params
my $filter = Genome::Model::Tools::Sam::VarFilter->create(
    use_version      => 'r963',
    min_read_depth   => 3,
    max_read_depth   => 100,
    gap_win_size     => 30,
    min_map_qual_snp => 25,
    gap_nearby_size  => 10,
    input_var_format => 'pileup',
    input_var_file   => $pileup_var_file,
    snv_out_file     => $pileup_out_files[0],
    indel_out_file   => $pileup_out_files[1],
    filtered_snv_out_file   => $pileup_out_files[2],
    filtered_indel_out_file => $pileup_out_files[3],
);

isa_ok($filter,'Genome::Model::Tools::Sam::VarFilter');
ok($filter->execute,'pileup varfilter executed ok');

$filter = Genome::Model::Tools::Sam::VarFilter->create(
    use_version      => 'r963',
    input_var_format => 'vcf',
    input_var_file   => $vcf_var_file,
    snv_out_file     => $vcf_out_files[0],
    indel_out_file   => $vcf_out_files[1],
    filtered_snv_out_file   => $vcf_out_files[2],
    filtered_indel_out_file => $vcf_out_files[3],
);

isa_ok($filter,'Genome::Model::Tools::Sam::VarFilter');
ok($filter->execute,'Vcf varfilter executed ok');

for my $format (@formats) {
    for my $i (0..$#out_files) {
        my $file_name = $out_files[$i];
        my $ori_file  = $root_dir.'/'.$file_name.'.'.$format;
        my $test_file = $tmp_dir.'/'.$file_name.'.'.$format;
        is(compare($ori_file, $test_file), 0, "$file_name.$format is generated as expected");
    }
}

done_testing();
