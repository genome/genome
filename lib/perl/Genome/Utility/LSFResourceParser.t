use strict;
use warnings FATAL => qw(all);
use Test::More;
use Data::Dumper;

BEGIN {
    use_ok('Genome::Utility::LSFResourceParser',
        'parse_lsf_params', 'construct_lsf_param_string')
}

sub parse_ok {
    my ($lsf_param_string, $lsf_params) = @_;
    my $parsed = parse_lsf_params($lsf_param_string);
    is_deeply($parsed, $lsf_params, "Parse: $lsf_param_string") or diag(Dumper({
        got => $parsed, expected => $lsf_params }));

    my $reconstructed_lsf_resource = construct_lsf_param_string($parsed);
    my $reparsed = parse_lsf_params($reconstructed_lsf_resource);
    is_deeply($reparsed, $lsf_params, "Unparse: $lsf_param_string") or
        diag(Dumper({got => $reparsed, expected => $lsf_params }));
}


parse_ok('', { 'options' => {}, 'rLimits' => {} });
parse_ok('rusage[mem=4000] span[hosts=1]', {
        'options' => {
            'resReq' => 'rusage[mem=4000] span[hosts=1]',
        },
        'rLimits' => {}
    });
parse_ok('rusage[tmp=100]', {
        'options' => {
            'resReq' => 'rusage[tmp=100]'
        },
        'rLimits' => {}
    });
parse_ok("-R 'span[hosts=1] rusage[mem=1000]' -n 4", {
        'options' => {
            'numProcessors' => '4',
            'resReq' => 'span[hosts=1] rusage[mem=1000]'
        },
        'rLimits' => {}
    });
parse_ok("-R 'span[hosts=1] rusage[mem=1000]' -n 4,6", {
        'options' => {
            'numProcessors' => '4',
            'maxNumProcessors' => '6',
            'resReq' => 'span[hosts=1] rusage[mem=1000]'
        },
        'rLimits' => {}
    });
parse_ok('select[tmp>1000] rusage[tmp=1000]', {
        'options' => {
            'resReq' => 'select[tmp>1000] rusage[tmp=1000]'
        },
        'rLimits' => {}
    });
parse_ok('select[mem>6000] rusage[mem=6000]', {
        'options' => {
            'resReq' => 'select[mem>6000] rusage[mem=6000]',
        },
        'rLimits' => {}
    });
parse_ok('rusage[mem=2000] select[mem > 2000] span[hosts=1]', {
        'options' => {
            'resReq' => 'rusage[mem=2000] select[mem > 2000] span[hosts=1]',
        },
        'rLimits' => {}
    });
parse_ok("-M 6000000 -R 'select[mem>16000] rusage[mem=16000]'", {

        'options' => {
            'resReq' => 'select[mem>16000] rusage[mem=16000]'
        },
        'rLimits' => {
            'RSS' => '6000000'
        }

    });
parse_ok(q{-R 'select[mem>16000] rusage[mem=16000]' -M 16000000}, {
        'options' => {
            'resReq' => 'select[mem>16000] rusage[mem=16000]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }
    });
parse_ok('select[mem>1024] rusage[mem=1024]', {
        'options' => {
            'resReq' => 'select[mem>1024] rusage[mem=1024]'
        },
        'rLimits' => {}
    });
parse_ok('rusage[mem=4000,tmp=1000] select[tmp>1000] span[hosts=1]', {
        'options' => {
            'resReq' => 'rusage[mem=4000,tmp=1000] select[tmp>1000] span[hosts=1]'
        },
        'rLimits' => {}
    });
parse_ok('rusage[mem=2000]', {

        'options' => {
            'resReq' => 'rusage[mem=2000]'
        },
        'rLimits' => {}
    });
parse_ok('-R \'select[mem>16000] rusage[mem=16000]\' -M 16000000 ', {
        'options' => {
            'resReq' => 'select[mem>16000] rusage[mem=16000]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }
    });
parse_ok('-M 16777216 rusage[mem=16384] select[mem > 16384] span[hosts=1]', {
        'options' => {
            'resReq' => '-M 16777216 rusage[mem=16384] select[mem > 16384] span[hosts=1]'
        },
        'rLimits' => {}

    });
parse_ok("-R 'select[mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]' -M 32000000", {
        'options' => {
            'resReq' => 'select[mem>32000 && gtmp>200] rusage[mem=32000:gtmp=200] span[hosts=1]'
        },
        'rLimits' => {
            'RSS' => '32000000'
        }
    });
parse_ok("-R 'select[mem>30000] rusage[mem=30000] span[hosts=1]' -M 30000000", {
        'options' => {
            'resReq' => 'select[mem>30000] rusage[mem=30000] span[hosts=1]'
        },
        'rLimits' => {
            'RSS' => '30000000'
        }
    });
parse_ok("-R 'select[gtmp>20 && mem>4000] span[hosts=1] rusage[gtmp=20,mem=4000]'", {
        'options' => {
            'resReq' => 'select[gtmp>20 && mem>4000] span[hosts=1] rusage[gtmp=20,mem=4000]'
        },
        'rLimits' => {}
    });
parse_ok("-M 8000000 -R 'select[mem>8000] rusage[mem=8000]'", {
        'options' => {
            'resReq' => 'select[mem>8000] rusage[mem=8000]'
        },
        'rLimits' => {
            'RSS' => '8000000'
        }
    });
parse_ok("-M 16000000 -R 'select[mem>16000] rusage[mem=16000]'", {
        'options' => {
            'resReq' => 'select[mem>16000] rusage[mem=16000]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }
    });
parse_ok("-M 12000000 -R 'select[mem>12000] rusage[mem=12000]'", {
        'options' => {
            'resReq' => 'select[mem>12000] rusage[mem=12000]'
        },
        'rLimits' => {
            'RSS' => '12000000'
        }
    });
parse_ok(q{-R 'select[mem>32000] rusage[mem=32000]' -M 32000000}, {
        'options' => {
            'resReq' => 'select[mem>32000] rusage[mem=32000]'
        },
        'rLimits' => {
            'RSS' => '32000000'
        }
    });
parse_ok(q{-R 'select[mem>12000] rusage[mem=12000]' -M 12000000}, {
        'options' => {
            'resReq' => 'select[mem>12000] rusage[mem=12000]'
        },
        'rLimits' => {
            'RSS' => '12000000'
        }
    });
parse_ok('select[mem>4096] rusage[mem=4096]', {

        'options' => {
            'resReq' => 'select[mem>4096] rusage[mem=4096]'
        },
        'rLimits' => {}
    });
parse_ok('rusage[tmp=2000] select[tmp>2000]', {

        'options' => {
            'resReq' => 'rusage[tmp=2000] select[tmp>2000]'
        },
        'rLimits' => {}

    });
parse_ok('rusage[mem=8000, tmp=2000] select[mem > 8000 && tmp > 2000] span[hosts=1] -M 8000000', {

        'options' => {
            'resReq' => 'rusage[mem=8000, tmp=2000] select[mem > 8000 && tmp > 2000] span[hosts=1] -M 8000000'
        },
        'rLimits' => {}

    });
parse_ok('rusage[mem=6000,tmp=10000] select[mem>6000 && tmp>10000] span[hosts=1]', {
        'options' => {
            'resReq' => 'rusage[mem=6000,tmp=10000] select[mem>6000 && tmp>10000] span[hosts=1]'
        },
        'rLimits' => {}
    });
parse_ok('-R \'select[mem>8192] rusage[mem=8192,tmp=100]\' -M 8192000 ', {

        'options' => {
            'resReq' => 'select[mem>8192] rusage[mem=8192,tmp=100]'
        },
        'rLimits' => {
            'RSS' => '8192000'
        }
    });
parse_ok('-R \'select[mem>32000] rusage[mem=32000]\' -M 32000000 ', {

        'options' => {
            'resReq' => 'select[mem>32000] rusage[mem=32000]'
        },
        'rLimits' => {
            'RSS' => '32000000'
        }
    });
parse_ok('-R "rusage[tmp=100]" -n 1', {
        'options' => {
            'numProcessors' => '1',
            'resReq' => 'rusage[tmp=100]'
        },
        'rLimits' => {}
    });
parse_ok('-R "rusage[tmp=100]" ', {
        'options' => {
            'resReq' => 'rusage[tmp=100]'
        },
        'rLimits' => {}

    });
parse_ok("-q short -R 'rusage[tmp=100]'", {
        'options' => {
            'queue' => 'short',
            'resReq' => 'rusage[tmp=100]'
        },
        'rLimits' => {}

    });
parse_ok("-q long rusage[tmp=100]", {
        'options' => {
            'resReq' => '-q long rusage[tmp=100]'
        },
        'rLimits' => {}
    });
parse_ok("-R 'span[hosts=1] rusage[mem=8000]' -M 8000000", {

        'options' => {
            'resReq' => 'span[hosts=1] rusage[mem=8000]'
        },
        'rLimits' => {
            'RSS' => '8000000'
        }

    });
parse_ok("-R 'select[tmp>2000] rusage[tmp=2000]'", {

        'options' => {
            'resReq' => 'select[tmp>2000] rusage[tmp=2000]'
        },
        'rLimits' => {}
    });
parse_ok("-R 'select[tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]' -M 16000000", {
        'options' => {
            'resReq' => 'select[tmp>1000 && mem>16000] span[hosts=1] rusage[tmp=1000:mem=16000]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }
    });
parse_ok("-R 'select[mem>=16000] rusage[mem=16000] span[hosts=1]' -M 16000000", {

        'options' => {
            'resReq' => 'select[mem>=16000] rusage[mem=16000] span[hosts=1]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }

    });
parse_ok("-R 'select[mem>32000] span[hosts=1] rusage[mem=32000]' -M 32000000 -n 2", {

        'options' => {
            'numProcessors' => '2',
            'resReq' => 'select[mem>32000] span[hosts=1] rusage[mem=32000]'
        },
        'rLimits' => {
            'RSS' => '32000000'
        }

    });
parse_ok("-R 'select[mem>32000 && tmp>50000] span[hosts=1] rusage[mem=32000,tmp=50000]' -M 32000000 -n 2", {
        'options' => {
            'numProcessors' => '2',
            'resReq' => 'select[mem>32000 && tmp>50000] span[hosts=1] rusage[mem=32000,tmp=50000]'
        },
        'rLimits' => {
            'RSS' => '32000000'
        }

    });
parse_ok("-R 'select[mem>16000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=16000]' -M 16000000", {

        'options' => {
            'resReq' => 'select[mem>16000 && tmp>150000] span[hosts=1] rusage[tmp=150000, mem=16000]'
        },
        'rLimits' => {
            'RSS' => '16000000'
        }
    });
parse_ok("-R 'select[mem>12000] rusage[mem=12000]' -M 12000000", {
        'options' => {
            'resReq' => 'select[mem>12000] rusage[mem=12000]'
        },
        'rLimits' => {
            'RSS' => '12000000'
        }
    });
parse_ok("-R 'select[mem=8192] rusage[mem=8192,tmp=1024]' -M 8192000", {
        'options' => {
            'resReq' => 'select[mem=8192] rusage[mem=8192,tmp=1024]'
        },
        'rLimits' => {
            'RSS' => '8192000'
        }
    });
parse_ok("-R 'select[gtmp>1] span[hosts=1] rusage[gtmp=1]'", {
        'options' => {
            'resReq' => 'select[gtmp>1] span[hosts=1] rusage[gtmp=1]'
        },
        'rLimits' => {}
    });
parse_ok("-R 'select[gtmp>1000] rusage[gtmp=1000] span[hosts=1]'", {

        'options' => {
            'resReq' => 'select[gtmp>1000] rusage[gtmp=1000] span[hosts=1]'
        },
        'rLimits' => {}

    });
parse_ok("-M 7000000 -R 'select[mem>7000 && tmp>10240] rusage[mem=7000]'", {
        'options' => {
            'resReq' => 'select[mem>7000 && tmp>10240] rusage[mem=7000]'
        },
        'rLimits' => {
            'RSS' => '7000000'
        }
    });
parse_ok("-M 5000000 -R 'select[mem>4500 && tmp>20000] rusage[mem=4500]'", {
        'options' => {
            'resReq' => 'select[mem>4500 && tmp>20000] rusage[mem=4500]'
        },
        'rLimits' => {
            'RSS' => '5000000'
        }
    });
parse_ok("-M 4000000 -R 'select[mem>4000] rusage[mem=4000]'", {
        'options' => {
            'resReq' => 'select[mem>4000] rusage[mem=4000]'
        },
        'rLimits' => {
            'RSS' => '4000000'
        }

    });
parse_ok("-M 25000000 -R 'select[mem>25000] rusage[mem=25000]'", {

        'options' => {
            'resReq' => 'select[mem>25000] rusage[mem=25000]'
        },
        'rLimits' => {
            'RSS' => '25000000'
        }
    });
parse_ok("-M 20000000 -R 'select[mem>20000] rusage[mem=20000]'", {
        'options' => {
            'resReq' => 'select[mem>20000] rusage[mem=20000]'
        },
        'rLimits' => {
            'RSS' => '20000000'
        }
    });
parse_ok("-M 200000 -n 4 -c 10 -R 'rusage[mem=200:gtmp=5]' -q short", {
        'options' => {
            'numProcessors' => '4',
            'queue' => 'short',
            'resReq' => 'rusage[mem=200:gtmp=5]'
        },
        'rLimits' => {
            'RSS' => '200000',
            'cpuTime' => '10'
        }
    });
parse_ok("-M 14000000 -R 'select[mem>14000] rusage[mem=14000]'", {
        'options' => {
            'resReq' => 'select[mem>14000] rusage[mem=14000]'
        },
        'rLimits' => {
            'RSS' => '14000000'
        }
    });

done_testing();
