#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 2;

my $subject = Genome::Sys->current_user();
my $view = Genome::Sys::User::View::Solr::Xml->create(subject => $subject);
ok($view->content, 'generated content');
ok(!($view->__errors__), 'view should not have errors') or diag explain $view->__errors__;
