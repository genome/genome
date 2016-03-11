#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use File::Basename qw();
use Genome::Test::Factory::Model::ReferenceSequence;

use Test::More tests => 4;
use Test::Deep qw(cmp_bag);

my $pkg = 'Genome::Model::Build::ReferenceSequence::Buckets';
use_ok($pkg);

my $reference = Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build;

my $seqdict_file = $reference->sequence_dictionary_path('sam');
my (undef, $seqdict_dir) = File::Basename::fileparse($seqdict_file);
Genome::Sys->create_directory($seqdict_dir);
Genome::Sys->write_file($seqdict_file, <<'EOFILE'
@HD	VN:1.0	SO:unsorted
@SQ	SN:1	LN:100	UR:file:///dev/null	AS:test1.0	M5:fakeuuid1	SP:test
@SQ	SN:2	LN:90	UR:file:///dev/null	AS:test1.0	M5:fakeuuid2	SP:test
EOFILE
);

my $result = $pkg->create(
    reference_sequence_build => $reference
);
isa_ok($result, $pkg, 'created buckets');

is($result->count, 2, 'count is set as expected');

my $bucket_list = $result->bucket_list;
cmp_bag($bucket_list, [[1],[2]], 'buckets created for each chromosome');
