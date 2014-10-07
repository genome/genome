#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require File::Spec;
require File::Temp;
use Test::Exception;
use Test::More;

my $class = 'Genome::Model::Tools::CgHub::Query';
use_ok($class) or die;

# Failures 
throws_ok(sub{ $class->execute(uuid => 'INVALID'); }, qr/\QFailed to find matching objects for 'INVALID' on CG Hub!\E/, 'create w/ invalid uuid fails');

# Success
my $uuid = '387c3f70-46e9-4669-80e3-694d450f2919';
my $tempdir = File::Temp::tempdir(CLEANUP => 1);
my $xml_file = File::Spec->catfile($tempdir, 'metadata.xml');
my $query = $class->create(
    uuid => $uuid,
    xml_file => $xml_file,
);
ok($query, 'create cg hub query cmd');
is($query->_build_command, "cgquery -o $xml_file analysis_id=$uuid", 'correct command');
is($query->uuid, $uuid, 'uuid');
ok($query->execute, 'execute cg hub query cmd');
ok(-s $query->xml_file, 'xml file downloaded'); # should be temp file

done_testing();
