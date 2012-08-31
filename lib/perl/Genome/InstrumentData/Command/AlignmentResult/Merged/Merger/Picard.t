#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";

test_default_max_jvm_heap_size();

done_testing();

sub test_default_max_jvm_heap_size {
    require Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard;
    require Genome::Sys;

    { # case 1: max_kb > safe_mem_limit_kb
        *Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::max_gb = sub { 6 };
        *Genome::Sys::mem_limit_kb = sub { 5*1_048_576 };
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::max_gb();
        ok($default_max_jvm_heap_size < $max_gb, '(case 1) default_max_jvm_heap_size is less than max_gb');
    }

    { # case 2: max_kb <= safe_mem_limit_kb
        *Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::max_gb = sub { 6 };
        *Genome::Sys::mem_limit_kb = sub { 6291456 };
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::max_gb();
        ok($default_max_jvm_heap_size < $max_gb, '(case 2) default_max_jvm_heap_size is less than max_gb');
    }

    { # case 3: no mem_limit_kb
        *Genome::Sys::mem_limit_kb = sub { 15728640 };
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::Merger::Picard::max_gb();
        ok($default_max_jvm_heap_size == $max_gb, '(case 3) default_max_jvm_heap_size equals max_gb');
    }
}

