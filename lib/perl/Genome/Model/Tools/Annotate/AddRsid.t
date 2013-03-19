use strict;
use warnings;

use above "Genome";
use File::Temp;
use Genome::Utility::Test qw(compare_ok);
use Test::More tests => 6;

my $class = 'Genome::Model::Tools::Annotate::AddRsid';
use_ok($class);

my $data_dir = Genome::Utility::Test->data_dir_ok($class);

my $tmp_file = Genome::Sys->create_temp_file_path();
my $cmd = Genome::Model::Tools::Annotate::AddRsid->create(
    anno_file   => "$data_dir/test.annotate.top2",
    vcf_file    => "$data_dir/test.vcf.gz",
    output_file => $tmp_file,
);

ok($cmd, 'created command object');
ok($cmd->execute, 'executed command');
ok(-s $tmp_file, 'output_file has size');
compare_ok($tmp_file, "$data_dir/test.out", 'output_file matched expected');
