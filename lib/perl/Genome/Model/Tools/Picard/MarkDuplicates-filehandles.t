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
    use_version => "1.34"
    );

my $fh = $cmd->get_max_filehandles_param;
like($fh, qr/^MAX_FILE_HANDLES=[0-9]+$/, "Version 1.34 sets MAX_FILE_HANDLES");

my $txt = "";
for my $v (1..33) {
    $cmd->use_version("1.$v");
    $txt .= $cmd->get_max_filehandles_param;
}
is($txt, "", "Versions < 1.34 do not set MAX_FILE_HANDLES");

done_testing();
