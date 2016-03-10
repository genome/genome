#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use File::Basename qw();
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::SoftwareResult::User;

use Test::More tests => 4;
use Test::Deep qw(cmp_bag);

my $pkg = 'Genome::Model::ReferenceSequence::Command::CreateBuckets';
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

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash;

my $cmd = $pkg->create(
    reference_sequence_build => $reference,
    result_users => $result_users,
);
isa_ok($cmd, $pkg, 'created command to create buckets');

ok($cmd->execute, 'executed command to created buckets');

isa_ok($cmd->output_result, $cmd->result_class, 'result created by command');

