#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above "Genome"; 
use Test::More;

use_ok('Genome::ProcessingProfile::View::Status::Xml') or die;

class Genome::Model::Tester { is => 'Genome::ModelDeprecated', };

my $pp = Genome::ProcessingProfile::Tester->create(
    name => 'Tester Test for Testing',
);
ok($pp, "created processing profile") or die;
my $model = Genome::Model->create(
    processing_profile => $pp,
    subject_name => 'human',
    subject_type => 'species_name',
);
ok($model, 'create model') or die;

my $view_obj = $pp->create_view(perspective => 'status', toolkit => 'xml'); 
ok($view_obj, "created a view") or die;
isa_ok($view_obj, 'Genome::ProcessingProfile::View::Status::Xml');

my $xml = $view_obj->_generate_content();
ok($xml, "view returns XML");

done_testing();
