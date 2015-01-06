#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::InstrumentData::Gatk::Base') or die;

class Genome::InstrumentData::Gatk::BaseTester {
    is => 'Genome::InstrumentData::Gatk::Base',
};

use Genome::Sys;
Sub::Install::reinstall_sub({
        code => sub { return; },
        into => 'Genome::Sys',
        as   => 'mem_limit_kb',
    });
is(Genome::InstrumentData::Gatk::BaseTester->max_memory_for_gmt_gatk, 4, 'default max memory for gmt gatk');
Sub::Install::reinstall_sub({
        code => sub { return 6000000; },
        into => 'Genome::Sys',
        as   => 'mem_limit_kb',
    });
is(Genome::InstrumentData::Gatk::BaseTester->max_memory_for_gmt_gatk, 5.7, 'max memory for gmt gatk');

my @instrument_data = (
    Genome::InstrumentData::Imported->__define__(id => -321),
);
my $alpha = Genome::InstrumentData::AlignmentResult::Merged->__define__();
$alpha->add_instrument_data_id(
    name => 'instrument_data_id',
    value_id => $instrument_data[0]->id,
    value_class_name => $instrument_data[0]->class,
);

my $bravo = Genome::InstrumentData::Gatk::BaseTester->__define__(bam_source => $alpha);
my $charlie = Genome::InstrumentData::Gatk::BaseTester->__define__(bam_source => $bravo);

is_deeply([$alpha->instrument_data], \@instrument_data, 'alpha instrument data');
is_deeply($bravo->bam_source, $alpha, 'bravo bam source is alpha');
is_deeply([$bravo->instrument_data], \@instrument_data, 'bravo instrument data');
is_deeply($charlie->bam_source, $bravo, 'charlie bam source is bravo');
is_deeply([$charlie->instrument_data], \@instrument_data, 'charlie instrument data');


done_testing();
