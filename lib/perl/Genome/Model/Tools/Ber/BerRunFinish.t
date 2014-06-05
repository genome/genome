#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

use Test::More tests => 3;

BEGIN {
    use_ok('Genome::Model::Tools::Ber::BerRunFinish');
    
}

my $tool_db = Genome::Model::Tools::Ber::BerRunFinish->create(
            'locus_tag'      => "PRABACTJOHNDFTBMGTEST",
		    'outdirpath'     => "/gscmnt/temp110/analysis/ktmp/BER_TEST/hmp/autoannotate.090715.old/out",
		    'sqlitedatafile' => "sqlite-PRABACTJOHNDFTBMGTEST-090417f.dat",
		    'sqliteoutfile'  => "sqlite-PRABACTJOHNDFTBMGTEST-090417f.out",
		    'acedb_version'  => "V3",
		    'amgap_path'     => Genome::Model::Tools::Hgmi->installation_path,
		    'assembly_version'     => "Version_1.0",
		    'pipe_version'   => "Version_1.0",
		    'project_type'   => "HGMI",
		    'org_dirname'     => "P_johnsonii",
		    'assembly_name'   => "Parabacteroides_johnsonii_PRABACTJOHNDFTBMGTEST_1.0.newb",
		    'sequence_set_id' => "221",
		);
isa_ok($tool_db,'Genome::Model::Tools::Ber::BerRunFinish');


SKIP: {
    skip "writes to annotation dir", 1;
    ok($tool_db->execute,'execute berrunfinish');
};

1;
