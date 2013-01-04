#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::Generator') or die;

class TestGenerator {
    is => 'Genome::Report::Generator',
};
sub TestGenerator::_add_to_report_xml{ return 1; };
sub TestGenerator::description { return 'Test report generator'; };

my $generator = TestGenerator->create();
ok($generator, 'create generator');

is($generator->name, 'Test Generator', 'name');
isa_ok($generator->generator, 'TestGenerator');
ok($generator->date, 'date');

done_testing();
