#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT}               = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

#$::RD_TRACE = 1;
#$::RD_HINT = 1;

use Test::More;

use above "Genome";

my $package = 'Genome::InstrumentData::Composite::Strategy';
use_ok($package)
  or die('test cannot continue');

my $strategy_fail = Genome::InstrumentData::Composite::Strategy->create(strategy =>
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4] v1'
);
isa_ok($strategy_fail, 'Genome::InstrumentData::Composite::Strategy', 'created strategy');
$strategy_fail->dump_status_messages(1);
ok(!$strategy_fail->execute, 'strategy parsing failed as expected');


_test_strategy(
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4] api v1',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
);

_test_strategy(
    'instrument_data '
    . 'aligned to contamination_ref using bwa 0.5.5 [-t 4] '
    . 'then filtered using dusting v1 '
    . 'api v1',
    {
        'action' => [
            {
                'params' => '',
                'parent' => {
                    'params'    => '-t 4',
                    'reference' => 'contamination_ref',
                    'version'   => '0.5.5',
                    'name'      => 'bwa',
                    'type'      => 'align'
                },
                'version' => 'v1',
                'name'    => 'dusting',
                'type'    => 'filter'
            }
        ],
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
);

_test_strategy(
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
        then (
                filtered using unaligned v1 then aligned to protein using rtg-mapx 1.2.3 then filtered using aligned v1
                and
                filtered using aligned v1
            )
            then aligned to virome_reference using rtg-map 1.2.3
     api v1',
    {
        'action' => [
            {
                'params' => '',
                'parent' => {
                    'params' => '',
                    'parent' => {
                        'params' => '',
                        'parent' => {
                            'params' => '',
                            'parent' => {
                                'params'    => '-t 4',
                                'reference' => 'contamination_ref',
                                'version'   => '0.5.5',
                                'name'      => 'bwa',
                                'type'      => 'align'
                            },
                            'version' => 'v1',
                            'name'    => 'unaligned',
                            'type'    => 'filter'
                        },
                        'reference' => 'protein',
                        'version'   => '1.2.3',
                        'name'      => 'rtg-mapx',
                        'type'      => 'align'
                    },
                    'version' => 'v1',
                    'name'    => 'aligned',
                    'type'    => 'filter'
                },
                'reference' => 'virome_reference',
                'version'   => '1.2.3',
                'name'      => 'rtg-map',
                'type'      => 'align'
            },
            {
                'params' => '',
                'parent' => {
                    'params' => '',
                    'parent' => {
                        'params'    => '-t 4',
                        'reference' => 'contamination_ref',
                        'version'   => '0.5.5',
                        'name'      => 'bwa',
                        'type'      => 'align'
                    },
                    'version' => 'v1',
                    'name'    => 'aligned',
                    'type'    => 'filter'
                },
                'reference' => 'virome_reference',
                'version'   => '1.2.3',
                'name'      => 'rtg-map',
                'type'      => 'align'
            }
        ],
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
);

_test_strategy(
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
        then (
                filtered using unaligned v1 then aligned to protein using rtg-mapx 1.2.3 then filtered using aligned v1
            )
            then aligned to virome_reference using rtg-map 1.2.3
     api v1',
    {
        'action' => [
        {
            'params' => '',
            'parent' => {
                'params' => '',
                'parent' => {
                    'params' => '',
                    'parent' => {
                        'params' => '',
                        'parent' => {
                            'params' => '-t 4',
                            'reference' => 'contamination_ref',
                            'version' => '0.5.5',
                            'name' => 'bwa',
                            'type' => 'align'
                        },
                        'version' => 'v1',
                        'name' => 'unaligned',
                        'type' => 'filter'
                    },
                    'reference' => 'protein',
                    'version' => '1.2.3',
                    'name' => 'rtg-mapx',
                    'type' => 'align'
                },
                'version' => 'v1',
                'name' => 'aligned',
                'type' => 'filter'
            },
            'reference' => 'virome_reference',
            'version' => '1.2.3',
            'name' => 'rtg-map',
            'type' => 'align'
        }
        ],
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
);


my $parent = {
    'params' => '-t 4',
    'reference' => 'contamination_ref',
    'version' => '0.5.5',
    'name' => 'bwa',
    'type' => 'align'
};
_test_strategy(
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
        then (
                filtered using unaligned v1 then aligned to protein using rtg-mapx 1.2.3 then filtered using aligned v1
                and
                filtered using aligned v1
                and
                filtered using dusting v1
            )
            then aligned to virome_reference using rtg-map 1.2.3
     api v1',
    {
        'action' => [
        {
            'params' => '',
            'parent' => {
                'params' => '',
                'parent' => {
                    'params' => '',
                    'parent' => {
                        'params' => '',
                        'parent' => $parent,
                        'version' => 'v1',
                        'name' => 'unaligned',
                        'type' => 'filter'
                    },
                    'reference' => 'protein',
                    'version' => '1.2.3',
                    'name' => 'rtg-mapx',
                    'type' => 'align'
                },
                'version' => 'v1',
                'name' => 'aligned',
                'type' => 'filter'
            },
            'reference' => 'virome_reference',
            'version' => '1.2.3',
            'name' => 'rtg-map',
            'type' => 'align'
        },
        {
            'params' => '',
            'parent' => {
                'params' => '',
                'parent' => $parent,
                'version' => 'v1',
                'name' => 'aligned',
                'type' => 'filter'
            },
            'reference' => 'virome_reference',
            'version' => '1.2.3',
            'name' => 'rtg-map',
            'type' => 'align'
        },
        {
            'params' => '',
            'parent' => {
                'params' => '',
                'parent' => $parent,
                'version' => 'v1',
                'name' => 'dusting',
                'type' => 'filter'
            },
            'reference' => 'virome_reference',
            'version' => '1.2.3',
            'name' => 'rtg-map',
            'type' => 'align'
        }
        ],
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
);

_test_strategy(
    'instrument_data 
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
     then merged using picard 1.29 then deduplicated using picard 1.29
     api v1',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'then' => {
            'params' => '',
            'then' => {
                'params' => '',
                'version' => '1.29',
                'name' => 'picard',
                'type' => 'deduplicate'
            },
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'data' => 'instrument_data',
        'api_version' => 'v1',
    },
    'parsed merge strategy as expected'
);

_test_strategy(
    'instrument_data
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
     then merged using picard 1.29 then deduplicated using picard 1.29
     then refined to variant_list using gatk-read-calibrator 0.01 [-et NO_ET]
     api v2',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'then' => {
            'params' => '',
            'then' => {
                'params' => '',
                'version' => '1.29',
                'name' => 'picard',
                'type' => 'deduplicate',
                then => {
                    params => '-et NO_ET',
                    version => '0.01',
                    name => 'gatk-read-calibrator',
                    type => 'refine',
                    known_sites => 'variant_list'
                }
            },
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'data' => 'instrument_data',
        'api_version' => 'v2',
    },
);

# Test new clip-overlap refiner
_test_strategy(
    'instrument_data
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
     then merged using picard 1.29 then deduplicated using picard 1.29
     then refined using clip-overlap 1.0.11
     api v2',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'then' => {
            'params' => '',
            'then' => {
                'params' => '',
                'version' => '1.29',
                'name' => 'picard',
                'type' => 'deduplicate',
                then => {
                    'params' => '',
                    'version' => '1.0.11',
                    'name' => 'clip-overlap',
                    'type' => 'refine',
                },
            },
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'data' => 'instrument_data',
        'api_version' => 'v2',
    },
);

# Test new and multiple refiners
_test_strategy(
    'instrument_data
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
     then merged using picard 1.29 then deduplicated using picard 1.29
     then refined to variant_list using gatk-read-calibrator 0.01 [-et NO_ET]
     then refined using clip-overlap 1.0.11
     api v2',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'then' => {
            'params' => '',
            'then' => {
                'params' => '',
                'version' => '1.29',
                'name' => 'picard',
                'type' => 'deduplicate',
                then => {
                    params => '-et NO_ET',
                    version => '0.01',
                    name => 'gatk-read-calibrator',
                    type => 'refine',
                    known_sites => 'variant_list',
                    then => {
                        'params' => '',
                        'version' => '1.0.11',
                        'name' => 'clip-overlap',
                        'type' => 'refine',
                    },
                }
            },
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'data' => 'instrument_data',
        'api_version' => 'v2',
    },
);

# Test new and multiple refiners in reverse order
_test_strategy(
    'instrument_data
     aligned to contamination_ref using bwa 0.5.5 [-t 4]
     then merged using picard 1.29 then deduplicated using picard 1.29
     then refined using clip-overlap 1.0.11
     then refined to variant_list using gatk-read-calibrator 0.01 [-et NO_ET]
     api v2',
    {
        'action' => [
            {
                'params'    => '-t 4',
                'reference' => 'contamination_ref',
                'version'   => '0.5.5',
                'name'      => 'bwa',
                'type'      => 'align'
            }
        ],
        'then' => {
            'params' => '',
            'then' => {
                'params' => '',
                'version' => '1.29',
                'name' => 'picard',
                'type' => 'deduplicate',
                then => {
                    'params' => '',
                    'version' => '1.0.11',
                    'name' => 'clip-overlap',
                    'type' => 'refine',
                    then => {
                        params => '-et NO_ET',
                        version => '0.01',
                        name => 'gatk-read-calibrator',
                        type => 'refine',
                        known_sites => 'variant_list',
                    },
                }
            },
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'data' => 'instrument_data',
        'api_version' => 'v2',
    },
    'parsed merge strategy as expected'
);

_test_strategy(
    'instrument_data aligned to reference using bwa 0.5.5 @festive-decoration [turkey] then merged using picard 1.29 api v3',
    {
        'api_version' => 'v3',
        'then' => {
            'params' => '',
            'version' => '1.29',
            'name' => 'picard',
            'type' => 'merge'
        },
        'action' => [{
            'params' => '',
            'decoration' => {
                'params' => 'turkey',
                'name' => 'festive-decoration'
            },
            'reference' => 'reference',
            'version' => '0.5.5',
            'name' => 'bwa',
            'type' => 'align'
        }],
       'data' => 'instrument_data'
    },
);

_test_strategy(
    'instrument_data both aligned to reference and merged using speedseq test_version api v3',
    {
        'api_version' => 'v3',
        'action' => [{
            'params' => '',
            'reference' => 'reference',
            'version' => 'test_version',
            'name' => 'speedseq',
            'type' => 'align_and_merge'
        }],
        'data' => 'instrument_data'
    }
);

_test_strategy(
    'instrument_data both aligned to reference and merged using speedseq test_version @festive-decoration [turkey] api v3',
    {
        'api_version' => 'v3',
        'action' => [{
            'params' => '',
            'decoration' => {
                'params' => 'turkey',
                'name' => 'festive-decoration'
            },
            'reference' => 'reference',
            'version' => 'test_version',
            'name' => 'speedseq',
            'type' => 'align_and_merge'
        }],
        'data' => 'instrument_data'
    }
);

done_testing();


sub _test_strategy {
    my $strategy_string = shift;
    my $expected_output = shift;

    my $strategy = $package->create(strategy => $strategy_string);
    isa_ok($strategy, $package, 'created strategy');
    ok($strategy->execute, 'execute strategy');
    is_deeply($strategy->tree, $expected_output, 'parsed strategy as expected')
        or diag Data::Dumper::Dumper($strategy->tree);

    return 1;
}
