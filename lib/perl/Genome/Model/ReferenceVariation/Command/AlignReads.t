#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ReferenceVariation;

use Sub::Override;

use Test::More tests => 4;

my $pkg = 'Genome::Model::ReferenceVariation::Command::AlignReads';
use_ok($pkg);

my $model = Genome::Test::Factory::Model::ReferenceVariation->setup_object;
for(1..3) {
    $model->add_instrument_data(
        Genome::Test::Factory::InstrumentData::Solexa->setup_object()
    );
}

my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

require Genome::SoftwareResult;

my $override = Sub::Override->new(
    'Genome::SoftwareResult::get_or_create',
    sub {
        package Genome::SoftwareResult;
        my $class = shift;
        return $class->SUPER::create;
    }
);

my $cmd = $pkg->create(build => $build);
isa_ok($cmd, $pkg, 'created command');

ok($cmd->execute, 'executed command');
like($cmd->status_message, qr(^Generated result), 'software result created');
