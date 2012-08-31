#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above 'Genome';

use_ok('Genome::InstrumentData::Solexa') or die;

my $pe2 = Genome::InstrumentData::Solexa->get(2862658358);
isa_ok($pe2,'Genome::InstrumentData::Solexa');
is($pe2->is_paired_end,1,'Paired End status found for lane');
is($pe2->calculate_alignment_estimated_kb_usage,300,'300kB disk needed for paired end instrument data');
my @fastq_files = @{$pe2->resolve_fastq_filenames};
is(scalar(@fastq_files),2,'got 2 fastq files for paired end instrument data');
# need to see if we can get the forward-only or reverse-only bases from the paird end inst data
is($pe2->total_bases_read('forward-only'),51200, 'forward only total_bases_read on paired end instrument data');
is($pe2->total_bases_read('reverse-only'),51200, 'reverse only total_bases_read on paired end instrument data');
is($pe2->total_bases_read('forward-only') + $pe2->total_bases_read('reverse-only'), $pe2->total_bases_read,
   'forward and reverse pairs add up to total bases');
is($pe2->read_count, 1024, 'read count');

done_testing();
exit;
