#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More;
use File::Compare;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 4;
}

use_ok('Genome::Model::Tools::Sam::Idxstats');

my @ori_map_chr_list = (1..22, 'X', 'Y', 'NT_113940');
my @full_chr_list = (@ori_map_chr_list, 'NT_113908', 'NT_113917');
my $sam_version = 'r963';

my $data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sam-Idxstats';
my $bam_file = $data_dir.'/test.bam';

my $tmp_dir  = Genome::Sys->create_temp_directory(Genome::Sys->username . "Idxstats_XXXXXX");
my $out_file = $tmp_dir . '/test.idxstats';

my $idxstat = Genome::Model::Tools::Sam::Idxstats->create(
    bam_file    => $bam_file,
    output_file => $out_file, 
    use_version => $sam_version,
);

isa_ok($idxstat,'Genome::Model::Tools::Sam::Idxstats');
my $result = $idxstat->execute;
ok($result, 'Samtools idxstats runs ok using sam version: ');

my $map_chr_list = $idxstat->map_ref_list($out_file);
my @final_map_chr_list; 

for my $chr_name (@full_chr_list) {
    push @final_map_chr_list, $chr_name if grep{$chr_name eq $_}@$map_chr_list;
}

is_deeply(\@ori_map_chr_list, \@final_map_chr_list, "Mapped chr list from idxstat using sam version $sam_version are created ok");

done_testing();
