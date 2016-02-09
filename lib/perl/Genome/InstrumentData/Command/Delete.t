#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use strict;
use warnings;

use Test::More tests => 4;

use above 'Genome';

use_ok('Genome::InstrumentData::Command') or die; # auto-gens delete command
ok(Genome::InstrumentData::Command::Delete->class, 'found delete command class') or die;

subtest 'delete imported instrument data' => sub{
    plan tests => 2;

    local $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;

    my $imported = Genome::InstrumentData::Imported->create(
        sequencing_platform => 'solexa',
        import_format => 'bam',
    );
    my $rv = Genome::InstrumentData::Command::Delete->_execute_with_shell_params_and_return_exit_code(
        '--instrument-data', $imported->id
    );
    is($rv, 0, 'successful return value when deleting imported instrument data');
    isa_ok($imported, 'UR::DeletedRef');

};

subtest 'do not delete solexa instrument data' => sub{
    plan tests => 2;

    local $ENV{UR_NO_REQUIRE_USER_VERIFY} = 1;

    my $solexa = Genome::InstrumentData::Solexa->create();
    my $rv = Genome::InstrumentData::Command::Delete->_execute_with_shell_params_and_return_exit_code(
        '--instrument-data', $solexa->id
    );
    is($rv, 1, 'failure return value when deleting solexa instrument data');
    isa_ok($solexa, 'Genome::InstrumentData::Solexa');

};

done_testing();
