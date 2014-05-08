#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 2;

require Genome::Env;

subtest 'set_default_values' => sub {
    plan tests => 4;

    no warnings 'once', 'redefine';

    local $ENV{GENOME_WITH_DEFAULT};
    local $ENV{GENOME_WITHOUT_DEFAULT};

    local *Genome::Env::allowed_modules = sub {
        return (
            'Genome::Env::GENOME_WITH_DEFAULT',
            'Genome::Env::GENOME_WITHOUT_DEFAULT',
        );
    };

    local *Genome::Env::GENOME_WITH_DEFAULT::default_value = sub {
        return 42;
    };

    ok(!defined $ENV{GENOME_WITHOUT_DEFAULT}, 'GENOME_WITHOUT_DEFAULT is not set');
    ok(!defined $ENV{GENOME_WITH_DEFAULT}, 'GENOME_WITH_DEFAULT is not set');

    Genome::Env::set_default_values();

    is($ENV{GENOME_WITH_DEFAULT}, 42, 'GENOME_WITH_DEFAULT got set');
    ok(!defined $ENV{GENOME_WITHOUT_DEFAULT}, 'GENOME_WITHOUT_DEFAULT is not set');
};

subtest 'check_genome_variables' => sub {
    plan tests => 2;

    no warnings 'once', 'redefine';

    local $ENV{GENOME_NOT_ALLOWED} = 1;

    do {
        local *STDERR;
        open(STDERR, '>', '/dev/null');

        local *Genome::Env::allowed_modules = sub {};
        ok(!Genome::Env::check_genome_variables(),
            'unrecognized variable causes check_genome_variables to return false');
    };

    my @allowed_modules = Genome::Env::allowed_modules();
    *Genome::Env::allowed_modules = sub {
        'Genome::Env::GENOME_NOT_ALLOWED',
        @allowed_modules,
    };
    ok(Genome::Env::check_genome_variables(),
        'recognized variable causes check_genome_variables to return true');
};
