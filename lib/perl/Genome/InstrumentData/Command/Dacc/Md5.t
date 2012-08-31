#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use Test::More;

use_ok('Genome::InstrumentData::Command::Dacc::Md5') or die;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Command-Dacc/SRS000000';

# BAM
my @bam_names = (qw/ BAM.1.bam /);
my @bam_files = map { $dir.'/'.$_ } @bam_names;
my $bam_md5 = Genome::InstrumentData::Command::Dacc::Md5->create(
    data_files => \@bam_files,
    format => 'bam',
    confirmed_md5_file => $dir.'/BAM.confirmed.md5',
);
ok($bam_md5, 'BAM: create');
$bam_md5->dump_status_messages(1);

is($bam_md5->_sra_id, 'SRS000000', 'BAM: SRA ID');
is($bam_md5->_directory, $dir, 'Directory');

my %expected_md5_files = $bam_md5->expected_md5_files;
#print Dumper(\%expected_md5_files);
is_deeply(
    \%expected_md5_files, 
    { map { $_.'.md5' => $_ } @bam_names },
    'BAM: expected md5 files',
);

my %dacc_md5 = $bam_md5->_load_dacc_md5;
#print Dumper(\%dacc_md5);
is_deeply(
    \%dacc_md5,
    { 'BAM.1.bam' => '10a3588b2babc2a6092bbe635b04169c' },
    'DACC md5',
);

my %confirmed_md5 = $bam_md5->_load_confirmed_md5;
#print Dumper(\%confirmed_md5);
is_deeply(
    \%dacc_md5,
    { 'BAM.1.bam' => '10a3588b2babc2a6092bbe635b04169c' },
    'Confirmed md5',
);
ok($bam_md5->execute, 'BAM: execute');

# SFF
my @sff_names = (qw/ SFF.1.sff SFF.2.sff /);
my @sff_files = map { $dir.'/'.$_ } @sff_names;
my $sff_md5 = Genome::InstrumentData::Command::Dacc::Md5->create(
    data_files => \@sff_files,
    format => 'sff',
    confirmed_md5_file => $dir.'/SFF.confirmed.md5',
);
ok($sff_md5, 'SFF: create');
$sff_md5->dump_status_messages(1);

is($sff_md5->_sra_id, 'SRS000000', 'SFF: SRA ID');
is($sff_md5->_directory, $dir, 'SFF: Directory');

%expected_md5_files = $sff_md5->expected_md5_files;
#print Dumper(\%expected_md5_files);
is_deeply(
    \%expected_md5_files, 
    { map { my $md5 = $_; $md5 =~ s/sff$/md5/; $md5 => $_ } @sff_names },
    'SFF: expected md5 files',
);

%dacc_md5 = $sff_md5->_load_dacc_md5;
#print Dumper(\%dacc_md5);
is_deeply(
    \%dacc_md5,
    {
        'SFF.2.sff' => 'd67d2d5e21961faedd47c67191b17ec8',
        'SFF.1.sff' => '2b4901370197c553c9e18592ec142b81'
    },
    'SFF: DACC md5',
);

%confirmed_md5 = $sff_md5->_load_confirmed_md5;
#print Dumper(\%confirmed_md5);
is_deeply(
    \%dacc_md5,
    {
        'SFF.2.sff' => 'd67d2d5e21961faedd47c67191b17ec8',
        'SFF.1.sff' => '2b4901370197c553c9e18592ec142b81'
    },
    'SFF: confirmed md5',
);
ok($sff_md5->execute, 'SFF: execute');

# FASTQ
my @fastq_names = (qw/ FASTQ.1.fastq.bz2 FASTQ.2.fastq.bz2 FASTQ.singleton.fastq.bz2 /);
my @fastq_files = map { $dir.'/'.$_ } @fastq_names;
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $fastq_md5 = Genome::InstrumentData::Command::Dacc::Md5->create(
    data_files => \@fastq_files,
    format => 'fastq',
    confirmed_md5_file => $tmpdir.'/confirmed.md5',
);
ok($fastq_md5, 'FASTQ: create');
$fastq_md5->dump_status_messages(1);

is($fastq_md5->_sra_id, 'SRS000000', 'FASTQ: SRA ID');
is($fastq_md5->_directory, $dir, 'FASTQ: Directory');

%expected_md5_files = $fastq_md5->expected_md5_files;
#print Dumper(\%expected_md5_files);
is_deeply(
    \%expected_md5_files, 
    { 'SRS000000.md5' => [ map { s/\.bz2//; $_ } @fastq_names ], },
    'FASTQ: expected md5 files',
);

%dacc_md5 = $fastq_md5->_load_dacc_md5;
#print Dumper(\%dacc_md5);
is_deeply(
    \%dacc_md5,
    {
        'FASTQ.1.fastq' => '6a5ffc295c8591d9f554c0f39676ea26',
        'FASTQ.2.fastq' => 'ecac97b4b7057ca2f7075701ec5bf533',
        'FASTQ.singleton.fastq' => '0b1239a96d938c0b7eb23b8448d63f22',
    },
    'FASTQ: DACC md5',
);

ok($fastq_md5->execute, 'FASTQ: execute'); # generate is tested here

%confirmed_md5 = $fastq_md5->_load_confirmed_md5;
#print Dumper(\%confirmed_md5);
is_deeply(
    \%dacc_md5,
    {
        'FASTQ.1.fastq' => '6a5ffc295c8591d9f554c0f39676ea26',
        'FASTQ.2.fastq' => 'ecac97b4b7057ca2f7075701ec5bf533',
        'FASTQ.singleton.fastq' => '0b1239a96d938c0b7eb23b8448d63f22',
    },
    'FASTQ: confirmed md5',
);

done_testing();
exit;

=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2010 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut

