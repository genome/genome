#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require Genome::Utility::Test;
use Test::Exception;
use Test::More;

use_ok('Genome::Model::Tools::CgHub::Query') or die;

my $data_dir = Genome::Utility::Test->data_dir_ok('Genome::Model::Tools::CgHub');
my $xml_file = File::Spec->join($data_dir, 'metadata.xml');
my $success_output_file = File::Spec->join($data_dir, 'query.out');
my $fail_output_file = File::Spec->join($data_dir, 'query.fail');

# So we don't actually send a request
my $run_command_cnt = 0;
sub Genome::Model::Tools::CgHub::Query::_run_command { $run_command_cnt++; return 1; };

# Failures 
my $query = Genome::Model::Tools::CgHub::Query->create(
    uuid => 'INVALID',
    output_file => $fail_output_file,
);
ok($query, 'create w/ invalid uuid');
throws_ok(sub{ $query->execute; }, qr/Ran CG Hub command, but it was determined that it was not successful!/, 'execute w/ invalid uuid fails');

# Success
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
$query = Genome::Model::Tools::CgHub::Query->create(
    uuid => $uuid,
    output_file => $success_output_file,
    xml_file => $xml_file,
);
ok($query, 'create cg hub query cmd');
is($query->_build_command, "cgquery 'analysis_id=$uuid' -o $xml_file", 'correct command');
ok($query->execute, 'execute cg hub query cmd');

is($run_command_cnt, 2, 'run command was invoked correctly');

done_testing();
