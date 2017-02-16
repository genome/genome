#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Test::More tests => 2;

use Genome::Test::Factory::Model::ReferenceAlignment qw();

my $model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
my $subject = Genome::ModelGroup->create(name => 'modelgroup view test', members => [$model]);
my $view = Genome::ModelGroup::View::Solr::Xml->create(subject => $subject);
ok($view->content, 'generated content');
ok(!($view->__errors__), 'view should not have errors') or diag explain $view->__errors__;
