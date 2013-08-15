#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::DeNovoAssembly::SxReadProcessor') or die;

class Genome::Model::Tools::Sx::Test::Default {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [ param => {}, ],
};
class Genome::Model::Tools::Sx::Test::InsertSize {
    is => 'Genome::Model::Tools::Sx::Base',
    has => [ min => {}, max => {}, ],
};

my $sample = Genome::Sample->create(name => '__TEST_SAMPLE__');
ok($sample, 'create test sample');
my $library = Genome::Library->create(name => '__TEST_LIBRARY__', sample => $sample);
ok($library, 'create library');
my @instrument_data;
my @inst_data_attrs = (
    {
        original_est_fragment_size => 250,
        is_paired_end => 1,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 3000,
        is_paired_end => 1,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 8000,
        is_paired_end => 1,
        read_count => 10000,
        read_length => 100,
    },
    {
        original_est_fragment_size => 100000,
        is_paired_end => 2,
        read_count => 10000,
        read_length => 777,
    },
);
my $i = 0;
for my $inst_data_attr ( @inst_data_attrs, ) {
    push @instrument_data, Genome::InstrumentData::Imported->create(
        library => $library,
        %$inst_data_attr,
    );
    $instrument_data[$i]->{sx_result_params} = {
        instrument_data_id => $instrument_data[$i]->id,
        output_file_count => ( $instrument_data[$i]->is_paired_end ? 2 : 1 ),
        output_file_type => 'sanger',
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    };
    $i++;
}
is(@instrument_data, @inst_data_attrs, 'create inst data');

diag('SUCCESS (NO READ PROCESSOR)');
my $sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create();
ok($sx_processor, 'create sx read processor');
ok($sx_processor->determine_processing(@instrument_data), 'determine processing');
my $no_read_processor_processing = { condition => 'DEFAULT', processor => '', };
is_deeply($sx_processor->default_processing, $no_read_processor_processing, 'default processing');
is_deeply($sx_processor->additional_processings, [], 'no additional processings');
my @expected_final_sx_result_params;
for my $instrument_data ( @instrument_data ) {
    my $processing = $sx_processor->determine_processing_for_instrument_data($instrument_data);
    is_deeply(
        $processing,
        $no_read_processor_processing,
        'correct processing for instrument data ',
    );
    my $sx_result_params = $sx_processor->sx_result_params_for_instrument_data($instrument_data),
    my %expected_sx_result_params = %{$instrument_data->{sx_result_params}};
    $expected_sx_result_params{read_processor} = $no_read_processor_processing->{processor};
    is_deeply(
        $sx_result_params,
        \%expected_sx_result_params,
        'correct sx result params for instrument data ',
    );
    push @expected_final_sx_result_params, \%expected_sx_result_params;
    ok(!$sx_processor->merged_sx_result_params_for_instrument_data($instrument_data), 'no merged sx result params');
}
is_deeply(
    [$sx_processor->final_sx_result_params],
    [sort {$a->{instrument_data_id} cmp $b->{instrument_data_id}} @expected_final_sx_result_params],
    'final sx result params',
);

diag('SUCCESS (OLD WAY)');
$sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'test default --param 1',
);
ok($sx_processor, 'create sx read processor');
ok($sx_processor->determine_processing(@instrument_data), 'determine processing');
my $old_way_processing = { condition => 'DEFAULT', processor => 'test default --param 1', };
is_deeply($sx_processor->default_processing, $old_way_processing, 'default processing');
is_deeply($sx_processor->additional_processings, [], 'no additional processings');
@expected_final_sx_result_params = ();
for my $instrument_data ( @instrument_data ) {
    my $processing = $sx_processor->determine_processing_for_instrument_data($instrument_data);
    is_deeply(
        $processing,
        $old_way_processing,
        'correct processing for instrument data ',
    );
    my $sx_result_params = $sx_processor->sx_result_params_for_instrument_data($instrument_data),
    my %expected_sx_result_params = %{$instrument_data->{sx_result_params}};
    $expected_sx_result_params{read_processor} = $old_way_processing->{processor};
    is_deeply(
        $sx_result_params,
        \%expected_sx_result_params,
        'correct sx result params for instrument data ',
    );
    push @expected_final_sx_result_params, \%expected_sx_result_params;
    ok(!$sx_processor->merged_sx_result_params_for_instrument_data($instrument_data), 'no merged sx result params');
}
is_deeply(
    [$sx_processor->final_sx_result_params],
    [sort {$a->{instrument_data_id} cmp $b->{instrument_data_id}} @expected_final_sx_result_params],
    'final sx result params',
);

diag('SUCCESS (NEW WAY DEFAULT ONLY)');
$sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFAULT (test default --param 1, coverage 10X)',
);
ok($sx_processor->determine_processing(@instrument_data), 'determine processing');
ok($sx_processor, 'create sx read processor');
my $new_way_processing = { condition => 'DEFAULT', processor => 'test default --param 1', coverage => 10, };
is_deeply($sx_processor->default_processing, $new_way_processing, 'default processing');
is_deeply($sx_processor->additional_processings, [], 'no additional processings');
my $merged_sx_result_params = $sx_processor->merged_sx_result_params_for_instrument_data($instrument_data[0]);
for my $instrument_data ( @instrument_data ) {
    diag($instrument_data->id);
    my $processing = $sx_processor->determine_processing_for_instrument_data($instrument_data);
    is_deeply(
        $processing,
        $new_way_processing,
        'correct processing for instrument data ',
    );
    my $sx_result_params = $sx_processor->sx_result_params_for_instrument_data($instrument_data),
    my %expected_sx_result_params = %{$instrument_data->{sx_result_params}};
    $expected_sx_result_params{read_processor} = $new_way_processing->{processor};
    is_deeply(
        $sx_result_params,
        \%expected_sx_result_params,
        'correct sx result params for instrument data ',
    );
    is_deeply(
        $sx_processor->merged_sx_result_params_for_instrument_data($instrument_data),
        $merged_sx_result_params,
        'merged sx result params',
    );
}
is_deeply(
    [$sx_processor->final_sx_result_params],
    [ $merged_sx_result_params ],
    'final sx result params',
);

diag('SUCESS (NEW WAY, FULL TEST)');
$sx_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFAULT (test default --param 1) original_est_fragment_size <= 2.5 * read_length (DEFAULT, coverage 30X) original_est_fragment_size > 1000 and original_est_fragment_size <= 6000 (test insert-size --min 1001 --max 6000 | test default --param 0, coverage 20X) original_est_fragment_size > 6000 and original_est_fragment_size <= 10000 (test insert-size --min 6001 --max 10000) read_length == 777 (test default --param 100)',
);
ok($sx_processor->determine_processing(@instrument_data), 'determine processing');
is_deeply($sx_processor->default_processing, $old_way_processing, 'default processing');
ok($sx_processor, 'create sx read processor');
ok($sx_processor->parser, 'parser');
my $default_processing = { condition => 'DEFAULT', processor => 'test default --param 1', };
my @expected_processings = (
    { condition => [qw/ original_est_fragment_size <= 2.5 * read_length /], processor => $default_processing->{processor}, coverage => 30, },
    { condition => [qw/ original_est_fragment_size > 1000 and original_est_fragment_size <= 6000 /], processor => 'test insert-size --min 1001 --max 6000 | test default --param 0', coverage => 20, },
    { condition => [qw/ original_est_fragment_size > 6000 and original_est_fragment_size <= 10000 /], processor => 'test insert-size --min 6001 --max 10000', },
    { condition => [qw/ read_length == 777 /], processor => 'test default --param 100', },
);
is_deeply($sx_processor->default_processing, $default_processing, 'default processor');
is_deeply($sx_processor->additional_processings, \@expected_processings, 'additional read processors');
@expected_final_sx_result_params = ();
for ( my $i = 0; $i < @instrument_data; $i++ ) {
    my $instrument_data = $instrument_data[$i];
    diag($instrument_data->id);
    my $processing = $sx_processor->determine_processing_for_instrument_data($instrument_data);
    is_deeply(
        $processing,
        $expected_processings[$i],
        'correct processing for instrument data ',
    );
    my $sx_result_params = $sx_processor->sx_result_params_for_instrument_data($instrument_data),
    my %expected_sx_result_params = %{$instrument_data->{sx_result_params}};
    $expected_sx_result_params{read_processor} = $expected_processings[$i]->{processor};
    is_deeply(
        $sx_result_params,
        \%expected_sx_result_params,
        'correct sx result params for instrument data ',
    );
    if ( $processing->{coverage} ) {
        my $merged_sx_result_params = $sx_processor->merged_sx_result_params_for_instrument_data($instrument_data);
        is_deeply(
            $merged_sx_result_params,
            {
                instrument_data_id => [ $instrument_data->id ],
                read_processor => $processing->{processor},
                coverage => $processing->{coverage},
                output_file_count => 2,
                output_file_type => 'sanger',
                test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
            },
            'merged sx result params',
        ); 
        push @expected_final_sx_result_params, $merged_sx_result_params;
    }
    else {
        ok(!$sx_processor->merged_sx_result_params_for_instrument_data($instrument_data), 'no merged sx result params');
        push @expected_final_sx_result_params, $sx_result_params;
    }
}
is_deeply(
    [$sx_processor->final_sx_result_params],
    [sort { (ref($a->{instrument_data_id})?$a->{instrument_data_id}->[0]:$a->{instrument_data_id}) cmp 
            (ref($b->{instrument_data_id})?$b->{instrument_data_id}->[0]:$b->{instrument_data_id})} @expected_final_sx_result_params],
    'final sx result params',
);

# FAILS
# invalid sx command
my $failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFAULT (test deefault --param 1)',
);
ok(!$failed_processor, 'failed to create sx read processor w/ invalid sx command');

# no default
$failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFULT (test default --param 1)',
);
ok(!$failed_processor, 'failed to create sx read processor w/ no default');

# more than one default
$failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFAULT (test default --param 1) DEFAULT (test default2 --param 1)',
);
ok(!$failed_processor, 'failed to create sx read processor w/ multiple defaults');

$failed_processor = Genome::Model::DeNovoAssembly::SxReadProcessor->create(
    processor => 'DEFAULT (test default --param 1)',
);
ok(!eval{$failed_processor->sx_result_params_for_instrument_data($instrument_data[0])}, 'failed to get sx result params before determining sx results for all instrument data');
ok(!eval{$failed_processor->merged_sx_result_params_for_instrument_data($instrument_data[0])}, 'failed to get merged sx result params before determining sx results for all instrument data');

done_testing();
