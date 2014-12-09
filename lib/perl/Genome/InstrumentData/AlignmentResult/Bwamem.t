#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use File::Temp;
use Path::Class;
use Test::More;

# This test was auto-generated because './InstrumentData/AlignmentResult/Bwamem.pm'
# had no '.t' file beside it.  Please remove this test if you believe it was
# created unnecessarily.  This is a bare minimum test that just compiles Perl
# and the UR class.
use_ok('Genome::InstrumentData::AlignmentResult::Bwamem');

test_input_fastq_removal();

done_testing();

sub test_input_fastq_removal {
    my ($tmpdir) = Path::Class::Dir->new( File::Temp::tempdir( CLEANUP => 1 ) );

    diag("Creating test file system in : $tmpdir");
    my ($root, @files) = create_test_file_system($tmpdir);

    my $rv = check_file_system_existence($root);
    is($rv, 1, "'$root' exists on file system");

    my $class = 'Genome::InstrumentData::AlignmentResult::Bwamem';
    ok($class->_remove_input_fastqs(@files), "removing fastqs");

    $rv = check_file_system_existence($root);
    is($rv, 0, "'$root' no longer exists on file system");
}

sub check_file_system_existence {
    my ($dir) = @_;

    unless (-e $dir) {
        diag("'$dir' path is missing");
        return 0;
    }

    return 1;
}

sub create_test_file_system {
    my $tmp = shift;

    my $root = $tmp->subdir('anonymous0');
    ok($root->mkpath, "Created root path: $root");

    my @files = ();
    for my $i (1..2) {
        my $name = join('_', 's', 'unknown', $i, 'sequence') . '.txt';
        my $f = $root->file($name);
        ok($f->touch, "Creating file: $f");
        push(@files, $f);
    }

    return $root, @files;
}
