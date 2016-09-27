#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use YAML;
use Test::More;
use File::Compare;
use Genome::Utility::Test;

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 10;
}

my $pkg = 'Genome::Model::Tools::Sam::Flagstat';
use_ok($pkg);

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, 'v4');
my $tmp_dir  = Genome::Sys->create_temp_directory(Genome::Sys->username . "Flagstat_XXXXXX");

my $bam_file      = $data_dir.'/test.bam';
my $ori_data_file = $data_dir.'/test.ori.data.yaml';
my $ori_data      = YAML::LoadFile($ori_data_file);

my $new_sam_version = 'r963';
my $old_sam_version = 'r599';

for my $sam_version ($new_sam_version, $old_sam_version) {
    my $out_file     = $tmp_dir . '/test.flagstat.' . $sam_version;
    my $ori_out_file = $data_dir .'/test.flagstat.' . $sam_version;

    my $flagstat = Genome::Model::Tools::Sam::Flagstat->create(
        bam_file    => $bam_file,
        output_file => $out_file, 
        use_version => $sam_version,  
    );

    isa_ok($flagstat,'Genome::Model::Tools::Sam::Flagstat');
    my $result = $flagstat->execute;
    ok($result, 'Flagstat tool runs ok using sam version: '.$sam_version);
    cmp_ok(compare($out_file, $ori_out_file), '==', 0, "Flagstat file by sam version $sam_version created ok.");
    my $data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($out_file);
    is_deeply($data, $ori_data, "Metrics parsed from flagstat using sam version $sam_version are created ok");
}

done_testing();
