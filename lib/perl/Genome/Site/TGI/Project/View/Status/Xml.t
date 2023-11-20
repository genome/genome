#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome"; 

use Test::More skip_all => 'requires connection to LIMS db';
#use Test::More tests => 6;
use Genome::Site::TGI::Project;

use_ok('Genome::Site::TGI::Project::View::Status::Xml') or die "test cannot continue...";

#2302675 is H1N1 Sequencing for Drug Resistance
my $subject = Genome::Site::TGI::Project->get(2302675);
ok($subject, "found expected project subject") or die "test cannot continue...";

my $view_obj = $subject->create_view(perspective => 'status', toolkit => 'xml'); 
ok($view_obj, "created a view") or die "test cannot continue...";
isa_ok($view_obj, 'Genome::Site::TGI::Project::View::Status::Xml');

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
