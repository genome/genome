#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Test::More tests => 6;

use_ok('Genome::FeatureList::View::Status::Xml') or die "test cannot continue...";

my $subject = Genome::FeatureList->get('linus2111.gsc.wustl.edu-7944-1287583590-1287583603-2542-10002');
ok($subject, "found expected feature-list subject") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'status', toolkit => 'xml');
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::FeatureList::View::Status::Xml');

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
