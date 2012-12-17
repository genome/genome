#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

use Test::More tests => 3;

BEGIN {
    use_ok('Genome::Model::Tools::Ber::BerRunBtabhmmtab');
    
}

my $tool_db = Genome::Model::Tools::Ber::BerRunBtabhmmtab->create(
                    'locus_tag'       => "PNI0002DFT",
		    'fastadirpath'    => "/gscmnt/temp110/analysis/ktmp/BER_TEST/hmp/autoannotate.090715.old/data/genomes/PNI0002DFT/fasta",
		    'berdirpath'      => "/gscmnt/temp110/analysis/ktmp/BER_TEST/hmp/autoannotate.090715.old/data/genomes/PNI0002DFT/ber",
		    'bsubfiledirpath' => "/gscmnt/temp110/analysis/ktmp/BER_TEST/hmp/autoannotate.090715.old/data/genomes/PNI0002DFT/bsubERRfiles",
		    'hmmdirpath'      => "/gscmnt/temp110/analysis/ktmp/BER_TEST/hmp/autoannotate.090715.old/data/genomes/PNI0002DFT/hmm",
		    'srcdirpath'      => "/gscmnt/temp110/info/annotation/ktmp/BER_TEST/hmp/autoannotate.090715.old/src",
		);
isa_ok($tool_db,'Genome::Model::Tools::Ber::BerRunBtabhmmtab');


SKIP: {
    skip "sends a couple thousand jobs to the blades.  needs a smaller data set for testing.", 1;
    ok($tool_db->execute,'execute berrunbtabhmmtab');
};

1;
