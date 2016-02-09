#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::More;

use above 'Genome';

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory';
class Thing::WithWorkingDirectory {
    is => 'Command::V2',
    roles => [ $class ],
};

my $tmp_dir = '/tmp';
my $thing = Thing::WithWorkingDirectory->create(working_directory => $tmp_dir);
ok($thing, 'created thing');

is(
    $thing->get_working_bam_path_with_new_extension('foo.bam', 'bar', 'baz'),
    File::Spec->join($tmp_dir, join('.', 'foo', 'bar', 'baz', 'bam')),
    'insert_extension_into_bam_path'
);

done_testing();
