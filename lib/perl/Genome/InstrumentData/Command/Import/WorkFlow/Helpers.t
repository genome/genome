#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Helpers';
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;

my $instrument_data = Genome::InstrumentData::Imported->__define__(
    original_data_path => '/dir/file.1.fastq,/dir/file.2.fastq.gz',
);
ok($instrument_data, 'define instrument data');
my $allocation = Genome::Disk::Allocation->__define__(
    owner => $instrument_data,
    mount_path => '/tmp',
    group_subdirectory => 'info',
    allocation_path => '100',
);
ok($allocation, 'define allocation');

my $helpers = $class->get;
ok($helpers, 'get helpers');
is_deeply(
    [ $helpers->local_source_files_for_instrument_data($instrument_data) ],
    [ map { $allocation->absolute_path.'/file.'.$_.'.fastq' } (1..2) ],
    'local source files for instrument data',
);

done_testing();
