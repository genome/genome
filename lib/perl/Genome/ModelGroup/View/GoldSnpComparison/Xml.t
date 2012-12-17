#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome"; 

use Test::More tests => 6;

use_ok('Genome::ModelGroup::View::GoldSnpComparison::Xml') or die "test cannot continue...";

my $subject = Genome::ModelGroup->get(242);
ok($subject, "defined a model group") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'gold-snp-comparison', toolkit => 'xml'); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::ModelGroup::View::GoldSnpComparison::Xml');

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

1;
