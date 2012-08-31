#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome'; 

use Test::More tests => 5;

# TODO: set up a fake build instead of using a real one
my $subject = Genome::Model::Build->get(106124235);

use_ok('Genome::Model::Build::View::DbSnpConcordance::Html') or die "test cannot continue...";

ok($subject, "found expected build subject") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'db-snp-concordance', toolkit => 'html'); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::Model::Build::View::DbSnpConcordance::Html');

my $html = $view_obj->_generate_content();
open(F, ">foo.html");
print F $html;
close(F);
ok($html, "view returns HTML") or die "test cannot continue...";

unlink "foo.html";
done_testing;
