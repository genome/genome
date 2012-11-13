#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

use_ok('Genome::Model::Tools::Sx::Rename') or die;

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $in_fastq = $dir.'/rename.in.fastq';
ok(-s $in_fastq, 'in fastq');
my $example_fastq = $dir.'/rename.example.fastq';
ok(-s $example_fastq, 'example fastq');

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Fail
ok( # no match_and_replace
    Genome::Model::Tools::Sx::Rename->create->__errors__,
    'errors when create w/o match_and_replace',
);
ok( # invalid match_and_replace
    Genome::Model::Tools::Sx::Rename->create(matches => [ 'qr/foo/g=bar' ])->__errors__,
    'errors when create w/ invalid match_and_replace',
);

# Ok
my $renamer = Genome::Model::Tools::Sx::Rename->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
    matches => [ 'qr|#.*/1$|=.b1', 'qr|#.*/2$|=.g1' ], # convert to gc read naming convention
    first_only => 1,
);
ok($renamer, 'create renamer');
isa_ok($renamer, 'Genome::Model::Tools::Sx::Rename');
ok($renamer->execute, 'execute renamer');
is(File::Compare::compare($example_fastq, $out_fastq), 0, "renamed as expected");

#print "$tmp_dir\n"; <STDIN>;
done_testing();
exit;

