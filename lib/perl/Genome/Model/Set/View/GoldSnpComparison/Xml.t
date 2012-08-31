#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome"; 
use Test::More;

if ($ENV{UR_RUN_LONG_TESTS}) {
    plan tests => 6;
}
else {
    plan skip_all => 'This test takes 20 minutes to run, run with UR_RUN_LONG_TESTS to enable';
}

use_ok('Genome::Model::Set::View::GoldSnpComparison::Xml') or die "test cannot continue...";

#apipe-test is our prefix for test models
my $subject = Genome::Model->define_set(name => { operator => 'LIKE', value => 'apipe-test%'});
ok($subject, "defined a model set") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'gold-snp-comparison', toolkit => 'xml'); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::Model::Set::View::GoldSnpComparison::Xml');

my $xml = $view_obj->_generate_content();
ok($xml, "view returns XML") or die "test cannot continue...";

SKIP: {
    skip "No Xml.t.expected in place.",1;
    my @diff =
        grep { $_ !~ /generated-at/ }
        grep { /\w/ }
        Genome::Sys->diff_file_vs_text(__FILE__ . '.expected',$xml);
    
    is("@diff","","XML has no differences from expected value");
}
