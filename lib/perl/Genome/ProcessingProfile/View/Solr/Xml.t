#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

use Genome::Test::Factory::ProcessingProfile::ReferenceAlignment qw();

my $subject = Genome::Test::Factory::ProcessingProfile::ReferenceAlignment->setup_object();
my $view = Genome::ProcessingProfile::View::Solr::Xml->create(subject => $subject);
ok($view->content, 'generated content');
ok(!($view->__errors__), 'view should not have errors') or diag explain $view->__errors__;
