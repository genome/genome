#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome"; 

use Test::More tests => 5;

use_ok('Genome::Model::Build::Set::View::Status::Xml') or die "test cannot continue...";

#The number of builds this returns may vary from time to time
my $subject = Genome::Model::Build::Set->get(status => { operator => 'LIKE', value => 'Running'});
ok($subject, "defined a build set") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'status', toolkit => 'xml'); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::Model::Build::Set::View::Status::Xml');

my $xml = $view_obj->_generate_content();
ok($xml, "view returns XML") or die "test cannot continue...";
