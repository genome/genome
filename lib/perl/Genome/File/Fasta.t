#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 10;
use Data::Dumper;

#test use first and quit if it doesn't work
use_ok('Genome::File::Fasta') or die;

#make a stupid small fai and reference to test
my $fake_fai = <<FAI;
1	5
2	10
3	4
FAI

my $fake_fasta = <<FASTA;
>1
ACGTA
>2
CGTAGCTGAC
>3
TTTT
FASTA

#start doing tests
#Create a temp dir for results
my $temp_dir = Genome::Sys->create_temp_directory();
ok($temp_dir, "created temp directory: $temp_dir");

#create the fasta and fai files

my $fasta_fh = Genome::Sys->open_file_for_writing("$temp_dir/test.fa");
print $fasta_fh $fake_fasta;
ok($fasta_fh->close, "wrote fasta");

my $fai_fh = Genome::Sys->open_file_for_writing("$temp_dir/test.fa.fai");
print $fai_fh $fake_fai;
ok($fai_fh->close, "wrote fai index");

my $file_obj = Genome::File::Fasta->create(id => "$temp_dir/test.fa");
ok($file_obj, "created Genome::File::Fasta object");

my $faidx_index = $file_obj->faidx_index;
is($faidx_index, "$temp_dir/test.fa.fai", "found fai file as expected");

my @expected_chromosome_names = (1,2,3);
is($file_obj->chromosome_name_list, @expected_chromosome_names, "got back chromosome names");

my @expected_chunks = (
    [[1,1,4]],
    [[1,5,5], [2,1,3]],
    [[2,4,7]],
    [[2,8,10], [3,1,1]],
    [[3,2,4]]
);
my @generated_chunks = $file_obj->divide_into_chunks(5);
is_deeply(\@generated_chunks, \@expected_chunks, "genome chunked as expected");

my @specific_chunk = $file_obj->divide_into_chunks(5,4);
is_deeply(\@specific_chunk, $expected_chunks[3], "specific chunk returned as expected");

eval {
    $file_obj->divide_into_chunks(5,6);
};
ok($@,"Bad specific chunk number throws as expected");
