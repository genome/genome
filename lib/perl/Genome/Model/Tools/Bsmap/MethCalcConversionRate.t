#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Utility::Test qw(compare_ok abort);
use Test::More tests => 8;

eval {
    # create a class instance
    my $class = 'Genome::Model::Tools::Bsmap::MethCalcConversionRate';
    use_ok($class);

    # check test data files
    my $data_dir = Genome::Utility::Test->data_dir($class);
    ok(-d $data_dir, "data_dir exists: $data_dir") or abort;

    # check inputs
    my $snvs_file = "$data_dir/snvs.hq";
    ok(-s $snvs_file, 'snvs file exists: snvs_file') or abort;

    my $output_file = "$data_dir/out_file";
    ok(-s $output_file, 'Output file exists: output_file') or abort;

    # create and execute command
    my $tmpdir = Genome::Sys->create_temp_directory();
    my $cmd = $class->create(
        snvs_file   => $snvs_file,
        output_file => "$tmpdir/out_file",
    );

    ok($cmd, 'created command') or abort;
    ok($cmd->execute, 'executed command') or abort;

    # check outputs
    ok(-s "$data_dir/out_file", 'output_file has size');
    my %compare_args = (
        replace => [
            [ qr(\Q$data_dir\E) => 'TEST_INPUTS_DIR' ],
        ],
    );
    compare_ok("$tmpdir/out_file", "$data_dir/out_file", 'output_file matched expected', %compare_args);
};

