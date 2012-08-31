#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More skip_all => 'This does not work yet';  #tests => 3;

use above 'Genome';

use_ok('Genome::Model::Build::ReferenceAlignment::View::ReferenceCoverage::Html');

my $subject = Genome::Model::Build::ReferenceAlignment->get(102576025);
ok($subject, 'found expected build subject');

my $view_obj = Genome::Model::Build::ReferenceAlignment::View::ReferenceCoverage::Html->create(subject_id => 102576025);
isa_ok($view_obj,'Genome::Model::Build::ReferenceAlignment::View::ReferenceCoverage::Html');
exit;
