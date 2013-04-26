#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::DeNovoAssembly::SxReadProcessor') or die;

my $sample = Genome::Sample->create(name => '__TEST_SAMPLE__');
ok($sample, 'create test sample');
my $library = Genome::Library->create(name => '__TEST_LIBRARY__', sample => $sample);
ok($library, 'create library');
my @instrument_data;
my @inst_data_attrs = (
    {
        original_est_fragment_size => 250,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 3000,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 8000,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 100000,
        read_count => 10000,
        read_length => 777,
    },
);
for my $inst_data_attr ( @inst_data_attrs ) {
    push @instrument_data, Genome::InstrumentData::Imported->create(
        library => $library,
        %$inst_data_attr,
    );
}
is(@instrument_data, @inst_data_attrs, 'create inst data');

diag('SUCCESS (OLD WAY)');
my $processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    instrument_data => \@instrument_data,
    read_processor => 'trim default --param 1',
);
ok($processor, 'failed to create sx read processor');
is_deeply($processor->_default_read_processing, { condition => 'DEFAULT', processor => 'trim default --param 1', }, 'got default processing');
for my $instrument_data ( @instrument_data ) {
    is_deeply( # all are default processing
        $processor->determine_processing_for_instrument_data($instrument_data),
        { condition => 'DEFAULT', processor => 'trim default --param 1', },
        'got correct processing for inst data',
    );
}

diag('SUCCESS (NEW WAY DEFAULT ONLY)');
$processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    instrument_data => \@instrument_data,
    read_processor => 'DEFAULT (trim default --param 1, coverage 10X)',
);
ok($processor, 'failed to create sx read processor');
is_deeply($processor->_default_read_processing, { condition => 'DEFAULT', processor => 'trim default --param 1', coverage => 10, }, 'got default processing');
for my $instrument_data ( @instrument_data ) {
    is_deeply( # all are default processing
        $processor->determine_processing_for_instrument_data($instrument_data),
        { condition => 'DEFAULT', processor => 'trim default --param 1', coverage => 10, },
        'got correct processing for inst data',
    );
}

diag('SUCESS (NEW WAY, FULL TEST)');
$processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    instrument_data => \@instrument_data,
    read_processor => 'DEFAULT (trim default --param 1) original_est_fragment_size <= 2.5 * read_length (DEFAULT, coverage 30X) original_est_fragment_size > 1000 and original_est_fragment_size <= 6000 (trim insert-size --min 1001 --max 6000 then filter by-length --length 50, coverage 20X) original_est_fragment_size > 6000 and original_est_fragment_size <= 10000 (trim insert-size --min 6001 --max 10000) read_length == 777 (filter --param 1 --qual 30)',
);
ok($processor, 'create sx read processor');
ok($processor->parser, 'got parser');
is_deeply($processor->_default_read_processing, { condition => 'DEFAULT', processor => 'trim default --param 1', }, 'got default processor');
my $processings = [
    { condition => [qw/ original_est_fragment_size <= 2.5 * read_length /], processor => 'DEFAULT', coverage => 30, },
    { condition => [qw/ original_est_fragment_size > 1000 and original_est_fragment_size <= 6000 /], processor => 'trim insert-size --min 1001 --max 6000 then filter by-length --length 50', coverage => 20, },
    { condition => [qw/ original_est_fragment_size > 6000 and original_est_fragment_size <= 10000 /], processor => 'trim insert-size --min 6001 --max 10000', },
    { condition => [qw/ read_length == 777 /], processor => 'filter --param 1 --qual 30', },
];
is_deeply($processor->_read_processings, $processings, 'got read processors');
for ( my $i = 0; $i < @instrument_data; $i++ ) {
    is_deeply(
        $processor->determine_processing_for_instrument_data($instrument_data[$i]),
        $processings->[$i],
        'got correct processing for inst data',
    );
}

# FAILS
# no default
my $failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    instrument_data => \@instrument_data,
    read_processor => 'DEFULT trim default --param 1',
);
ok(!$failed_processor, 'failed to create sx read processor');

# more than one default
$failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    instrument_data => \@instrument_data,
    read_processor => 'DEFAULT (trim default --param 1) DEFAULT (trim default2 --param 1)',
);
ok(!$failed_processor, 'failed to create sx read processor');



done_testing();

