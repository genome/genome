#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Sub::Install;
use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::Query') or die;
use_ok('Genome::Model::Tools::CgHub::Test') or die;

# So we don't actually send a request
Genome::Model::Tools::CgHub::Test->overload_lwp_user_agent_request;

# Failures 
throws_ok(sub{ Genome::Model::Tools::CgHub::Query->execute(query => 'INVALID'); }, qr/Failed to execute query! INVALID/, 'execute w/ invalid query fails');

# Success
my $analysis_id = '387c3f70-46e9-4669-80e3-694d450f2919';
my $query = Genome::Model::Tools::CgHub::Query->execute(
    query => 'analysis_id='.$analysis_id,
);
ok($query->result, 'execute cg hub query cmd');
like($query->metadata_xml, qr/^\<\?xml/, 'metadata xml set');

done_testing();
