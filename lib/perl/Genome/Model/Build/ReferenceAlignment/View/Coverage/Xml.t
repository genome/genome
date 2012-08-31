#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml');

my $subject = Genome::Model::Build::ReferenceAlignment->get(102576025);
ok($subject, 'found expected build subject');

my $view_obj = Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml->create(
    subject_id => 102576025,
    aspects => ['status'],
);
isa_ok($view_obj,'Genome::Model::Build::ReferenceAlignment::View::Coverage::Xml');
my $xml = $view_obj->_generate_content;
ok($xml,'got xml content');
exit;
