#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Sub::Install;
use Test::Exception;
use Test::MockObject;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::Query') or die;

# So we don't actually send a request
my $request_cnt = 0;
my $response = Test::MockObject->new();
$response->mock('content', sub{ return "<XML>\n"; });
Sub::Install::reinstall_sub({
        code => sub{
            $request_cnt++;
            if ( ${$_[1]->uri} =~ /INVALID/ ) {
                $response->set_false('is_success');
            }
            else {
                $response->set_true('is_success');
            }
            return $response;
        },
        into => 'LWP::UserAgent',
        as => 'request',
    });

# Failures 
throws_ok(sub{ Genome::Model::Tools::CgHub::Query->execute(query => 'INVALID'); }, qr/Failed to execute query! INVALID/, 'execute w/ invalid query fails');

# Success
my $analysis_id = '387c3f70-46e9-4669-80e3-694d450f2919';
my $query = Genome::Model::Tools::CgHub::Query->execute(
    query => 'analysis_id='.$analysis_id,
);
ok($query->result, 'execute cg hub query cmd');
is($request_cnt, 2, 'LWP::UserAgent:request called 2X');

done_testing();
