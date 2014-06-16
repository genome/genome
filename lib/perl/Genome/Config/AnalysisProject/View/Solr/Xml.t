#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;
use above 'Genome';

use_ok('Genome::Config::AnalysisProject::View::Solr::Xml') or die 'failed to use module';

my $ap = Genome::Config::AnalysisProject->create(
    name => 'Test Project for ' . __FILE__
);
ok($ap, 'created a test Analysis Project');

my $view = $ap->create_view(perspective => 'solr', toolkit => 'xml');
isa_ok($view, 'Genome::Config::AnalysisProject::View::Solr::Xml', 'created view');

my $xml = $view->_generate_content();
ok($xml, 'view returns XML');


