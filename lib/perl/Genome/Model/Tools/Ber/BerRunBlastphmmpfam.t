#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

use Test::More tests => 3;

BEGIN {
    use_ok('Genome::Model::Tools::Ber::BerRunBlastphmmpfam');
    
}

my $tool_db = Genome::Model::Tools::Ber::BerRunBlastphmmpfam->create(
                    'locus_tag'       => "PNI0002DFT",
		    'proteinnamefile' => "PNI0002DFT.proteinname.fof",
		    'fastadirpath'    => "/gscmnt/temp110/analysis/wnash/BER_TEST/PNI0002DFT/fasta",
		    'berdirpath'      => "/gscmnt/temp110/analysis/wnash/BER_TEST/PNI0002DFT/ber",
		    'bsubfiledirpath' => "/gscmnt/temp110/analysis/wnash/BER_TEST/PNI0002DFT/bsubERRfiles",
		    'hmmdirpath'      => "/gscmnt/temp110/analysis/wnash/BER_TEST/PNI0002DFT/hmm",
		    'blpqueryfile'    => "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/data/panda/AllGroup/AllGroup.niaa",
		    'hmmdatabase'     => "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate/data/ALL_LIB.20081108.HMM",

		);
isa_ok($tool_db,'Genome::Model::Tools::Ber::BerRunBlastphmmpfam');


SKIP: {
    skip "sends a couple thousand jobs to the blades.  needs a smaller data set for testing.", 1;
    ok($tool_db->execute,'execute berrunblastphmmpfam');
};
