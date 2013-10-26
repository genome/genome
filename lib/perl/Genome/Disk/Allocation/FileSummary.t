#! /usr/bin/env genome-perl
use strict;
use warnings;

$ENV{UR_DBI_NO_COMMIT} = 1;

use Test::More;
use above "Genome";
use File::Basename qw(basename dirname);
use File::Spec;
use lib File::Spec->join(dirname(__FILE__), 't-lib');
use GenomeDiskAllocationCommon qw(create_test_volumes);

use_ok('Genome::Disk::Allocation') or die;
use_ok('Genome::Disk::Volume') or die;

my @volumes = create_test_volumes(2);

my $allocation = Genome::Disk::Allocation->create(
    disk_group_name => $volumes[0]->disk_group_names,
    allocation_path => 'testing123',
    owner_class_name => 'UR::Value',
    owner_id => 'foo',
    kilobytes_requested => 1,
);

my $tempdir = Genome::Sys->create_temp_directory();

my $real_symlink_destination = File::Spec->join($tempdir,'blah');
my $broken_symlink_destination = File::Spec->join($tempdir, 'broke');
my $directory_symlink_destination = File::Spec->join($tempdir, 'this_is_a_directory');
my $empty_file = 'empty_file';
my $non_empty_file = 'non_empty_file';
my $real_symlink = 'this_is_a_symlink';
my $broken_symlink = 'this_symlink_is_broken';
my $directory_symlink = 'symlink_to_a_dir';
my $test_directory = 'test_directory';

my $empty_file_path = File::Spec->join($allocation->absolute_path, $empty_file);
my $non_empty_file_path = File::Spec->join($allocation->absolute_path, $non_empty_file);
my $real_symlink_path = File::Spec->join($allocation->absolute_path, $real_symlink);
my $broken_symlink_path = File::Spec->join($allocation->absolute_path, $broken_symlink);
my $directory_symlink_path = File::Spec->join($allocation->absolute_path, $directory_symlink);
my $test_directory_path = File::Spec->join($allocation->absolute_path, $test_directory);

Genome::Sys->write_file($empty_file_path, '');
Genome::Sys->write_file($non_empty_file_path, 'turkey');
Genome::Sys->write_file($real_symlink_destination, 'sandwich');
Genome::Sys->create_directory($directory_symlink_destination);
Genome::Sys->create_directory($directory_symlink_destination);
Genome::Sys->create_directory($test_directory_path);
symlink($real_symlink_destination, $real_symlink_path);
symlink($broken_symlink_destination, $broken_symlink_path);
symlink($directory_symlink_destination, $directory_symlink_path);


#is live symlink
assert_values($allocation, $real_symlink, 1, 1, 1, $real_symlink_destination);

#is broken symlink
assert_values($allocation, $broken_symlink, 1, 0, 0, $broken_symlink_destination);

#is symlink to directory
assert_values($allocation, $directory_symlink, 1, 0, 0, $directory_symlink_destination);

#empty file
assert_values($allocation, $empty_file, 0, 1, 1, undef);

#non empty file
assert_values($allocation, $non_empty_file, 0, 1, 1, undef);

#die when given directory
eval {
    Genome::Disk::Allocation::FileSummary->create(
        allocation_id => $allocation->id,
        file => $test_directory,
    )
};

my $error = $@;
like($error, qr(is a directory), 'it should die if given a directory');

$allocation->delete();

done_testing();

sub assert_values {
    my $allocation = shift;
    my $file = shift,
    my $symlink_val = shift;
    my $digest_val = shift;
    my $size_val = shift;
    my $destination_val = shift;

    my $summary = Genome::Disk::Allocation::FileSummary->create(
        allocation_id => $allocation->id,
        file => $file,
    );

    is($summary->is_symlink, $symlink_val, 'it correctly determined if it was looking at a symlink');
    if ($digest_val) {
        ok(defined($summary->digest), 'it set an md5 digest value');
    } else {
        is($summary->digest, undef, 'it did not set an md5 digest value');
    }
    if ($size_val) {
        ok(defined($summary->size_in_bytes), 'it set a value for size')
    } else {
        is($summary->size_in_bytes, undef, 'it did not set a value for size');
    }
    is($summary->destination, $destination_val, 'destination pointed where it was supposed to');
}
