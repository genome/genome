#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

require File::Temp;
require Genome::Utility::Test;
use Test::Exception;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::Helpers';
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

my $helpers = $class->get;
ok($helpers, 'get helpers');

# source files functions
my @source_files = (
    $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Command-Import-Basic/fastq-1.fq.gz',
    $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-Command-Import-Basic/fastq-2.fastq',
);

# source file retrieval
throws_ok(sub{ $helpers->source_file_retrieval_method(); }, qr/No source file to get retrieval method!/, 'source_file_retrieval_method fails w/o source file');
my %source_files_and_retrieval_methods = (
    'https://cghub.ucsc.edu/UUID' => 'cg hub',
    'file://some/url.edu/file' => 'remote url',
    'http://some.url.edu/file' => 'remote url',
    '/some/local/file' => 'local disk',
);
for my $sf ( keys %source_files_and_retrieval_methods ) {
    is(
        $helpers->source_file_retrieval_method($sf),
        $source_files_and_retrieval_methods{$sf},
        "source file $sf retrieval method is '$source_files_and_retrieval_methods{$sf}'",
    );
}

throws_ok(
    sub{ $helpers->source_files_retrieval_method(); },
    qr/No source files to get retrieval method!/,
    'source_file_retrieval_methods fails w/o source file',
);
throws_ok(
    sub{ $helpers->source_files_retrieval_method(keys %source_files_and_retrieval_methods); },
    qr/Mixed file retrieval methods for source files! /,
    'source_files_retrieval_method fails w/ many retrieval methods',
);
my @source_files_from_remote_url = grep { $source_files_and_retrieval_methods{$_} eq 'remote url' } keys %source_files_and_retrieval_methods;
is(@source_files_from_remote_url, 2, '2 files from reote url');
is(
    $helpers->source_files_retrieval_method(@source_files_from_remote_url),
    'remote url',
    "source_file_retrieval_methods is 'remote url'",
);

# source file format
ok(!eval{$helpers->source_file_format()}, 'format for no source file fails');
ok(!$helpers->source_file_format('source.duh'), 'no format for unknown source file');
is($helpers->error_message, 'Unrecognized source file format! source.duh', 'correct error');
is($helpers->source_file_format($source_files[0]), 'fastq', 'source file 1 format');
is($helpers->source_file_format($source_files[1]), 'fastq', 'source file 2 format');
is($helpers->source_file_format('source_file.fastq.tgz'), 'fastq', 'format for tgz source file is fastq');
is($helpers->source_file_format('source_file.fastq.tar.gz'), 'fastq', 'format for tar.gz source file is fastq');
is($helpers->source_file_format('source.bam'), 'bam', 'format for bam source file is bam');
is($helpers->source_file_format('source.sra'), 'sra', 'format for sra source file is sra');
is($helpers->source_file_format('source.fasta'), 'fasta', 'format for fasta source file is fasta');
is($helpers->source_file_format('source.fa'), 'fasta', 'format for fa source file is fasta');
is($helpers->source_file_format('source.fna'), 'fasta', 'format for fna source file is fasta');

# is_source_file_archived
throws_ok(sub {$helpers->is_source_file_archived; }, qr/No source file to determined if archived!/, 'is_source_file_archived failed w/o source file');
ok($helpers->is_source_file_archived('file.tar.gz'), 'file.tar.gz is archived');
ok($helpers->is_source_file_archived('file.tar'), 'file.tar is archived');
ok($helpers->is_source_file_archived('file.tgz'), 'file.tgz is archived');
ok(!$helpers->is_source_file_archived('filetar.gz'), 'filetar.gz is not archived');

ok(!eval{$helpers->size_of_source_file;}, 'failed to get size for source file w/o source file');
ok(!eval{$helpers->size_of_remote_file;}, 'failed to get size for remote file w/o remote file');

ok(!eval{$helpers->kilobytes_required_for_processing_of_source_files;}, 'failed to get kilobytes needed for processing w/o source files');
is($helpers->kilobytes_required_for_processing_of_source_files(@source_files), 746, 'kilobytes needed for processing source files');

# headers
ok(!eval{$helpers->load_headers_from_bam;}, 'failed to get headers w/o bam');
my $input_bam = $test_dir.'/bam-rg-multi/v1/input.rg-multi.bam';
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
            'ID:2883581798	PL:illumina	PU:2883581798.	LB:TEST-patient1-somval_normal1-extlibs	PI:165	DS:paired end	DT:2012-12-17T13:15:46-0600	SM:TEST-patient1-somval_normal1	CN:WUGSC',
            'ID:2883581799	PL:illumina	PU:2883581799.	LB:TEST-patient1-somval_normal1-extlibs	PI:165	DS:paired end	DT:2012-12-17T13:15:46-0600	SM:TEST-patient1-somval_normal1	CN:WUGSC',
        ],
    },
    'headers',
);

ok(!eval{$helpers->read_groups_from_headers;}, 'failed to get read groups from headers w/o headers');
my $read_groups_from_headers = $helpers->read_groups_from_headers($headers->{'@RG'});
is_deeply(
    $read_groups_from_headers, 
    {
        2883581797 => 'CN:WUGSC	DS:paired end	DT:2012-12-17T13:15:46-0600	LB:TEST-patient1-somval_normal1-extlibs	PI:165	PL:illumina	PU:2883581797.	SM:TEST-patient1-somval_normal1',
        2883581798 => 'CN:WUGSC	DS:paired end	DT:2012-12-17T13:15:46-0600	LB:TEST-patient1-somval_normal1-extlibs	PI:165	PL:illumina	PU:2883581798.	SM:TEST-patient1-somval_normal1',
        2883581799 => 'CN:WUGSC	DS:paired end	DT:2012-12-17T13:15:46-0600	LB:TEST-patient1-somval_normal1-extlibs	PI:165	PL:illumina	PU:2883581799.	SM:TEST-patient1-somval_normal1',
    },
    'read groups from headers',
);

ok(!eval{$helpers->headers_to_string;}, 'failed headers to string w/o headers');
my $headers_string = $helpers->headers_to_string($headers);
ok($headers_string, 'headers to string');# cannot compare for some reason

ok(!eval{$helpers->load_read_groups_from_bam;}, 'failed to load read groups from bam w/o bam');
is_deeply($helpers->load_read_groups_from_bam($test_dir.'/sra/v1/input.sra.bam'), [], 'non read groups in bam');
is_deeply($helpers->load_read_groups_from_bam($input_bam), [qw/ 2883581797 2883581798 2883581799 /], 'load read groups from bam');

# verify tmp disk
ok($helpers->verify_adequate_disk_space_is_available_for_source_files(working_directory => '/tmp', source_files => \@source_files), 'verify adequate disk space is available for source files');

# flagstat
my $bam_basename = 'input.bam';
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $bam_path = $tmp_dir.'/'.$bam_basename;
Genome::Sys->create_symlink($test_dir.'/bam/v1/'.$bam_basename, $bam_path);
ok(-s $bam_path, 'linked bam path');
my $run_flagstat = $helpers->load_or_run_flagstat($bam_path); # runs
ok($run_flagstat, 'run flagstat');
my $load_flagstat = $helpers->load_or_run_flagstat($bam_path); # loads
is_deeply($load_flagstat, $run_flagstat, 'load flagstat');
ok($helpers->validate_bam($bam_path), 'validate bam');

# md5
throws_ok { $helpers->md5_path_for } qr/^No path given to get md5 path!/, 'failed md5_path_for undef';
is($helpers->md5_path_for($bam_path), $bam_path.'.md5', 'md5_path_for');
throws_ok { $helpers->original_md5_path_for } qr/^No path given to get original md5 path!/, 'failed original_data_path_md5 undef';
is($helpers->original_md5_path_for($bam_path), $bam_path.'.md5-orig', 'original_md5_path_for');

my $run_md5 = $helpers->load_or_run_md5($bam_path); # runs
ok($run_md5, 'run md5');
my $load_md5 = $helpers->load_or_run_md5($bam_path); # loads
is_deeply($load_md5, $run_md5, 'load md5');

# previously imported
my @md5s = map { $_ x 32 } (qw/ a b c /);
my @instrument_data = map {
    Genome::InstrumentData::Imported->__define__(
        id => $_ - 11,
    )
} (0..$#md5s);
is(@instrument_data, 3, '__define__ instdata');
ok(!$helpers->were_original_path_md5s_previously_imported(md5s => \@md5s), 'as expected, no inst data found for md5s');
my @mdr_attrs = map { 
    Genome::InstrumentDataAttribute->__define__(
    instrument_data_id => $instrument_data[$_]->id,
    attribute_label => 'original_data_path_md5',
    attribute_value => $md5s[$_],
    nomenclature => 'WUGC',
) } (0..$#md5s);
is(@mdr_attrs, 3, '__define__ md5 attrs');
my $ds_attr = Genome::InstrumentDataAttribute->__define__(
    instrument_data_id => $instrument_data[2]->id,
    attribute_label => 'downsample_ratio',
    attribute_value => 0.25,
    nomenclature => 'WUGC',
);
ok($ds_attr, '__define__ downsample attr for instdata 3');

## no downsample ratio
ok(# a & b w/o downsample_ratio should be found
    $helpers->were_original_path_md5s_previously_imported(md5s => \@md5s),
    'inst data found for md5s "a" & "b" w/o downsample_ratio',
);
is($helpers->error_message, 'Instrument data was previously imported! Found existing instrument data: -10, -11', 'correct error message');
ok(# c w/o downsample_ratio should not be found
    !$helpers->were_original_path_md5s_previously_imported(md5s => [$md5s[2]]),
    'as expected, no inst data found for "c" md5',
);

## w/ downsample ratio
ok(# a, b & c w/ downsample_ratio 0.33 should not be found
    !$helpers->were_original_path_md5s_previously_imported(md5s => \@md5s, downsample_ratio => 0.33),
    'as expected, no inst data found for "a", "b" & "c" md5 w/ downsample_ratio of .33',
);
ok(# a & b w/ downsample_ratio 0.25 should not be found
    !$helpers->were_original_path_md5s_previously_imported(md5s => [$md5s[0..1]], downsample_ratio => 0.25),
    'as expected, no inst data found for "a" and "b" md5 w/ downsample_ratio of .25',
);
ok(# c w/ downsample_ratio of 0.25 should be found
    $helpers->were_original_path_md5s_previously_imported(md5s => \@md5s, downsample_ratio => 0.25),
    'instdata found for md5 "c" and downsample_ratio 0.25',
);
is($helpers->error_message, 'Instrument data was previously downsampled by a ratio of 0.25 and imported! Found existing instrument data: -9', 'correct error message');

# properties
my $properties = $helpers->key_value_pairs_to_hash(qw/ sequencing_platform=solexa lane=2 flow_cell_id=XXXXXX /);
is_deeply(
    $properties,
    { sequencing_platform => 'solexa', lane => 2, flow_cell_id => 'XXXXXX', },
    'key value piars to hash',
);
$properties = $helpers->key_value_pairs_to_hash(qw/ sequencing_platform=solexa lane=2 lane=3 flow_cell_id=XXXXXX /);
ok(!$properties, 'failed as expected to convert key value pairsr into hash with duplicate label');
is($helpers->error_message, "Multiple values for instrument data property! lane => 2, 3", 'correct error');
$properties = $helpers->key_value_pairs_to_hash(qw/ sequencing_platform=solexa lane= flow_cell_id=XXXXXX /);
is($helpers->error_message, 'Failed to parse with instrument data property label/value! lane=', 'correct error');

# rm source files
ok(!eval{$helpers->remove_path_and_auxiliary_files();}, 'failed to remove source paths and md5s w/o source paths');
Genome::Sys->create_symlink($test_dir.'/bam/v1/'.$bam_basename.'.md5-orig', $bam_path.'.md5-orig');
Genome::Sys->create_symlink($test_dir.'/bam/v1/'.$bam_basename.'.md5-orig', $bam_path.'.random');
ok($helpers->remove_paths_and_auxiliary_files($bam_path), 'remove source paths and md5s w/o source paths');
ok(!glob($bam_path.'*'), 'removed path and auxillary files');

# bam path
throws_ok( sub{ $helpers->insert_extension_into_bam_path(); }, qr/^No bam path given to insert extension to bam path!/, 'insert_extension_into_bam_path fails w/o bam_path');
throws_ok( sub{ $helpers->insert_extension_into_bam_path('in.bam'); }, qr/^No extension given to insert extension to bam path!/, 'insert_extension_into_bam_path fails w/o extension');
throws_ok( sub{ $helpers->insert_extension_into_bam_path('bam', 'sorted'); }, qr/^Failed to insert extension into bam path! Bam path does not end in .bam! bam/, 'insert_extension_into_bam_path fails w/ invalid bam');
is($helpers->insert_extension_into_bam_path('in.bam', 'sorted'), 'in.sorted.bam', 'insert_extension_into_bam_path');

# validators
throws_ok(sub{ $helpers->is_downsample_ratio_invalid(); }, qr/No downsample ratio to check!/, 'is_downsample_ratio_invalid fails w/o downsample_ratio');

my @errors = $helpers->is_downsample_ratio_invalid('NA');
ok(@errors, 'errors for downsample_ratio of NA');
is($errors[0]->desc, 'Invalid number! NA', 'correct error desc for downsample_ratio of NA');

@errors = $helpers->is_downsample_ratio_invalid('0');
ok(@errors, 'errors for downsample_ratio of 0');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 0', 'correct error desc for downsample_ratio of 0');

@errors = $helpers->is_downsample_ratio_invalid('1');
ok(@errors, 'errors for downsample_ratio of 1');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 1', 'correct error desc for downsample_ratio of 1');

# is_bam_paired_end
throws_ok(sub{ $helpers->is_bam_paired_end(); }, qr/No bam path given to is_bam_paired_end!/, 'is_bam_paired_end fails w/o bam');
throws_ok(sub{ $helpers->is_bam_paired_end('does_not_exist'); }, qr/Bam path given to is_bam_paired_end does not exist!/, 'is_bam_paired_end fails w/ non existing bam');
my $data_dir2 = File::Spec->join($test_dir, 'bam-rg-multi', 'v4');
my $bam_path2 = File::Spec->join($data_dir2, '2883581797.paired.bam');
my $is_paired_end = $helpers->is_bam_paired_end($bam_path2);
is($is_paired_end, 1, "bam $bam_path2 is paired end");
$bam_path2 = File::Spec->join($data_dir2, '2883581797.singleton.bam');
$is_paired_end = $helpers->is_bam_paired_end($bam_path2);
is($is_paired_end, 0, "bam $bam_path2 is not paired end");

# update bam
my $instdata = Genome::InstrumentData::Imported->__define__(
    id => -1111,
);
ok($instdata, '__define__ instdata');
Sub::Install::reinstall_sub({
        code => sub{ return $bam_path2; },
        into => "Genome::InstrumentData::Imported",
        as   => 'bam_path',
    });
$instdata->add_attribute( # test removal
    attribute_label => 'read_length',
    attribute_value => -1,
    nomenclature => 'WUGC',
);
ok($ds_attr, '__define__ downsample attr for instdata 3');
throws_ok(sub{ $helpers->update_bam_metrics_for_instrument_data(); }, qr/No instrument data given to update bam for instrument data!/, 'update_bam_metrics_for_instrument_data fails w/o instdata');
ok($helpers->update_bam_metrics_for_instrument_data($instdata, $bam_path2), 'update_bam_metrics_for_instrument_data');
my %expected_attrs = (
    bam_path => $bam_path2,
    read_count => 94,
    read_length => 100,
    is_paired_end => 0,
);
for my $attr (qw/ bam_path is_paired_end read_count read_length /) {
    is($instdata->attributes(attribute_label => $attr)->attribute_value, $expected_attrs{$attr}, "$attr set")
}

done_testing();
