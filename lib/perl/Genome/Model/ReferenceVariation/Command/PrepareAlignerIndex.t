#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::ReferenceVariation;

use Sub::Override;

use Test::More tests => 4;

my $pkg = 'Genome::Model::ReferenceVariation::Command::PrepareAlignerIndex';
use_ok($pkg);

my $model = Genome::Test::Factory::Model::ReferenceVariation->setup_object;
my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

require Genome::Model::Build::ReferenceSequence::IndexBase;
my $index_override = Sub::Override->new(
    'Genome::Model::Build::ReferenceSequence::IndexBase::get_or_create',
    sub {
        my $class = shift;
        return $class->UR::Object::create;
    }
);

my $cmd = $pkg->create(build => $build);
isa_ok($cmd, $pkg, 'created command');

ok($cmd->execute, 'executed command');
like($cmd->status_message, qr(^Generated result), 'software result created');
