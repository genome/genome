#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More;

my $class = "Genome::Model::Tools::Picard::MarkDuplicates";
use_ok($class);


my $cmd = $class->create(
    input_file => '/tmp/foo',
    output_file => '/tmp/bar',
    metrics_file => '/tmp/baz',
    temp_directory => '/tmp',
    use_version => "1.36"
    );

my $fh = $cmd->build_cmdline_string;
like($fh, qr/MAX_FILE_HANDLES=[0-9]+$/, "Version > 1.34 sets MAX_FILE_HANDLES");

$cmd->use_version("1.31");
my $txt .= $cmd->build_cmdline_string;
unlike($txt, qr/MAX_FILE_HANDLES=[0-9]+$/, "Version < 1.34 does not set MAX_FILE_HANDLES");

done_testing();
