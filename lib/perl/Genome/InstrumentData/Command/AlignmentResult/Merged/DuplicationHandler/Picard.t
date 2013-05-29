#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use Sub::Install;
use above "Genome";

test_default_max_jvm_heap_size();

done_testing();

sub test_default_max_jvm_heap_size {
    require Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard;
    require Genome::Sys;

    { # case 1: max_kb > safe_mem_limit_kb
        Sub::Install::reinstall_sub({
            into => 'Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard',
            as => 'max_gb',
            code => sub { 6 }
        });
        Sub::Install::reinstall_sub({
            into => 'Genome::Sys',
            as => 'mem_limit_kb',
            code => sub { 5*1_048_576 },
        });
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::max_gb();
        ok($default_max_jvm_heap_size < $max_gb, '(case 1) default_max_jvm_heap_size is less than max_gb');
    }

    { # case 2: max_kb <= safe_mem_limit_kb
        Sub::Install::reinstall_sub({
            into => 'Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard',
            as => 'max_gb',
            code =>  sub { 6 },
        });
        Sub::Install::reinstall_sub({
            into => 'Genome::Sys',
            as => 'mem_limit_kb',
            code => sub { 6291456 },
        });
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::max_gb();
        ok($default_max_jvm_heap_size < $max_gb, '(case 2) default_max_jvm_heap_size is less than max_gb');
    }

    { # case 3: no mem_limit_kb
        Sub::Install::reinstall_sub({
            into => 'Genome::Sys',
            as => 'mem_limit_kb',
            code => sub { 15728640 },
        });
        my $default_max_jvm_heap_size = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::default_max_jvm_heap_size();
        my $max_gb = Genome::InstrumentData::Command::AlignmentResult::Merged::DuplicationHandler::Picard::max_gb();
        ok($default_max_jvm_heap_size == $max_gb, '(case 3) default_max_jvm_heap_size equals max_gb');
    }
}
