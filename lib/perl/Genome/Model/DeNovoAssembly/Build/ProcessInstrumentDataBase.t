#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Genome::Model::Build::DeNovoAssembly;
use Test::MockObject;
use Test::More;

my $class = 'Genome::Model::DeNovoAssembly::Build::ProcessInstrumentDataBase';
use_ok($class) or die;

my $lsf_resource = $class->lsf_resource;
ok($lsf_resource, 'lsf resource');

my $pp = Test::MockObject->new;
$pp->mock('read_processor', sub{ return 'trim far --version 1.7 --nr-threads 4';});
my $build = Test::MockObject->new;
$build->mock(
    'resolve_resource_requirements_for_processing_instrument_data',
    sub{ return Genome::Model::Build::DeNovoAssembly::resolve_resource_requirements_for_processing_instrument_data($build, $class); }
);
$build->mock('processing_profile', sub{ return $pp; });

my $resource_req =  $class->resolve_resource_requirements(
    {build => $build},
); 
is($resource_req, $lsf_resource.' -n 4', 'resource requirements');

done_testing();
