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

my $original_directory = File::Spec->join(File::Spec->rootdir, 'original');
my $working_directory = File::Spec->join(File::Spec->rootdir, 'working-directory');
my $thing = Thing::WithWorkingDirectory->create(working_directory => $working_directory);
ok($thing, 'created thing');

is(
    $thing->get_working_path_for_file_path( File::Spec->join($original_directory, 'tmp', 'foo.fastq') ),
    File::Spec->join($working_directory, 'foo.fastq'),
    'get_working_path_for_file_path'
);
is(
    $thing->get_working_bam_path_with_new_extension( File::Spec->join($original_directory, 'foo.bam'), 'bar', 'baz'),
    File::Spec->join($working_directory, join('.', 'foo', 'bar', 'baz', 'bam')),
    'insert_extension_into_bam_path'
);

done_testing();
