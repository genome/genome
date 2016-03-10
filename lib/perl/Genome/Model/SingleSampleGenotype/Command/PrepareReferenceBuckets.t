#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use File::Basename qw();
use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::SingleSampleGenotype;

use Test::More tests => 4;

my $pkg = 'Genome::Model::SingleSampleGenotype::Command::PrepareReferenceBuckets';
use_ok($pkg);

my $model = Genome::Test::Factory::Model::SingleSampleGenotype->setup_object;
my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

my $reference = $build->reference_sequence_build;
my $seqdict_file = $reference->sequence_dictionary_path('sam');
my (undef, $seqdict_dir) = File::Basename::fileparse($seqdict_file);
Genome::Sys->create_directory($seqdict_dir);
Genome::Sys->write_file($seqdict_file, <<'EOFILE'
@HD	VN:1.0	SO:unsorted
@SQ	SN:1	LN:100	UR:file:///dev/null	AS:test1.0	M5:fakeuuid1	SP:test
@SQ	SN:2	LN:90	UR:file:///dev/null	AS:test1.0	M5:fakeuuid1	SP:test
EOFILE
);

my $cmd = $pkg->create(
    build => $build,
);
isa_ok($cmd, $pkg, 'created command to prepare buckets');

ok($cmd->execute, 'executed command to prepare buckets');
ok($reference->buckets, 'buckets were created');

