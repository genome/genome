#! /gsc/bin/perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Helpers';
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

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

# source files functions
my @source_files = (
    $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Command-Import-Basic/fastq-1.txt.gz',
    $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Command-Import-Basic/fastq-2.fastq',
);
is($helpers->source_file_format($source_files[0]), 'fastq', 'source file 1 format');
is($helpers->source_file_format($source_files[1]), 'fastq', 'source file 2 format');

ok(!eval{$helpers->size_of_source_file;}, 'failed to get size for source file w/o source file');
ok(!eval{$helpers->size_of_remote_file;}, 'failed to get size for remote file w/o remote file');

ok(!eval{$helpers->kilobytes_required_for_processing_of_source_files;}, 'failed to get kilobytes needed for processing w/o source files');
is($helpers->kilobytes_required_for_processing_of_source_files(@source_files), 738, 'kilobytes needed for processing source files');

# headers
ok(!eval{$helpers->load_headers_from_bam;}, 'failed to get headers w/o bam');
my $input_bam = $test_dir.'/input.two-read-groups.bam';
ok(-s $input_bam, 'input bam');
my $headers = $helpers->load_headers_from_bam($input_bam);
is_deeply(
    $headers,
    {
        '@SQ' => [
            'SN:MT	LN:16299	UR:ftp://ftp.ncbi.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh37/special_requests/GRCh37-lite.fa.gz	AS:GRCh37-lite	M5:c68f52674c9fb33aef52dcf399755519	SP:Mus musculus'
        ],
        '@PG' => [
            'ID:2883581797	VN:0.5.9	CL:bwa aln -t4 -q 5; bwa sampe  -a 435.616825 ',
            'ID:2883581798	VN:0.5.9	CL:bwa aln -t4 -q 5; bwa sampe  -a 435.616825 ',
            'ID:GATK PrintReads	VN:f6dac2d	CL:readGroup=null platform=null number=-1 downsample_coverage=1.0 sample_file=[] sample_name=[] simplify=false no_pg_tag=false'
        ],
        '@HD' => [ 'VN:1.4	GO:none	SO:coordinate' ],
        '@RG' => [
            'ID:2883581797	PL:illumina	PU:2883581797.	LB:TEST-patient1-somval_normal1-extlibs	PI:165	DS:paired end	DT:2012-12-17T13:15:46-0600	SM:TEST-patient1-somval_normal1	CN:WUGSC',
            'ID:2883581798	PL:illumina	PU:2883581798.	LB:TEST-patient1-somval_normal1-extlibs	PI:165	DS:paired end	DT:2012-12-17T13:15:46-0600	SM:TEST-patient1-somval_normal1	CN:WUGSC'
        ]
    },
    'headers',
);

ok(!eval{$helpers->read_groups_from_headers;}, 'failed to get read groups from headers w/o headers');
my $read_groups_from_headers = $helpers->read_groups_from_headers($headers->{'@RG'});
is_deeply(
    $read_groups_from_headers, 
    {
        2883581797 => 'CN:WUGSC	DS:paired end	DT:2012-12-17T13:15:46-0600	LB:TEST-patient1-somval_normal1-extlibs	PI:165	PL:illumina	PU:2883581797.	SM:TEST-patient1-somval_normal1',
        2883581798 => 'CN:WUGSC	DS:paired end	DT:2012-12-17T13:15:46-0600	LB:TEST-patient1-somval_normal1-extlibs	PI:165	PL:illumina	PU:2883581798.	SM:TEST-patient1-somval_normal1'
    },
    'read groups from headers',
);

ok(!eval{$helpers->headers_to_string;}, 'failed headers to string w/o headers');
my $headers_string = $helpers->headers_to_string($headers);
ok($headers_string, 'headers to string');# cannot compare for some reason

ok(!eval{$helpers->load_read_groups_from_bam;}, 'failed to load read groups from bam w/o bam');
is_deeply([$helpers->load_read_groups_from_bam($input_bam)], [qw/ 2883581797 2883581798 /], 'load read groups from bam');

# verify tmp disk
ok($helpers->verify_adequate_disk_space_is_available_for_source_files(tmp_dir => '/tmp', source_files => \@source_files), 'verify adequate disk space is available for source files');

done_testing();
