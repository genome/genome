#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

use Test::More tests => 3;

BEGIN {
    use_ok('Genome::Model::Tools::Ber::BerRunAnnoSqlite');
    
}

my $tool_db = Genome::Model::Tools::Ber::BerRunAnnoSqlite->create(
                    'locus_tag'  => "PNI0002DFT",
		    'srcdirpath' => "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/src",
		    'outdirpath' => "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/out",
		);
isa_ok($tool_db,'Genome::Model::Tools::Ber::BerRunAnnoSqlite');

SKIP: {
    skip "writes to annotation dir",1; # should check if in gscana
    ok($tool_db->execute,'execute berrunannosqlite');
};
