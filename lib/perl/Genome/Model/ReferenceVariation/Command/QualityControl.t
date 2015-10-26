#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Model::ReferenceVariation;

use Sub::Override;

use Test::More tests => 5;

my $pkg = 'Genome::Model::ReferenceVariation::Command::QualityControl';
use_ok($pkg);

my $model = Genome::Test::Factory::Model::ReferenceVariation->setup_object;
for(1..3) {
    $model->add_instrument_data(
        Genome::Test::Factory::InstrumentData::Solexa->setup_object()
    );
}

my $build = Genome::Test::Factory::Build->setup_object(model_id => $model->id);

my $speedseq_result = Genome::Test::Factory::InstrumentData::AlignmentResult::Merged::Speedseq->setup_object(
    instrument_data => [$build->instrument_data],
    reference_build => $build->reference_sequence_build,
    map {; $_ => $model->$_ } (qw(aligner_name aligner_version aligner_params)),
);
$speedseq_result->add_user(user => $build, label => 'merged_alignment_result');
is($build->merged_alignment_result, $speedseq_result, 'alignment result added to build');

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
