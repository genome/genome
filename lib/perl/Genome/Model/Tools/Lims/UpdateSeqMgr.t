#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More; #skip_all => 'test data not in place yet....';

use_ok( 'Genome::Model::Tools::Lims::UpdateSeqMgr' );

#TODO - fix and move to test suite
#my $project_dir = '/gscmnt/936/info/jschindl/pbtestout';
#Genome::Model::Tools::Lims::UpdateSeqMgr->execute(project_dir => $project_dir);

done_testing();

exit;

