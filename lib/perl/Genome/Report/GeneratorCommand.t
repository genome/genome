#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Report::GeneratorCommand');
my %params = (
    print_xml => 1,
    #print_datasets => 1,
    #datasets => 'rows',
    #email => Genome::Config->user_email,
);
my $cmd = Genome::Report::GeneratorCommand->create(%params);

my $report = $cmd->_generate_report_and_execute_functions(
    name => 'Rows n Stuff',
    description => 'Testing the generator command',
    headers => [qw/ column1 column2 colum3 /],
    rows => [ [qw/ row1.1 row1.2 row1.3 /], [qw/ row2.1 row2.2 row2.3 /] ],
);
ok($report, 'Generated report');
isa_ok($report, 'Genome::Report');

done_testing();
exit;

