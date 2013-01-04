#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Compare;
use Test::More skip_all => "Does not play nice with the test harness";
#use Test::More tests => 7;

BEGIN {
        use_ok('Genome::Model::Tools::WuBlast::Parse');
}

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-WuBlast';

my @output_files = map{$test_dir."/blast$_.blast"}(1..3);
my @parse_files  = map{$test_dir."/parse$_.out"}(1..3);
map{unlink $_ if -e $_}@parse_files;

my $parse1 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $output_files[0],
    parse_outfile => $parse_files[0],
);

ok($parse1->execute, 'Parse1 finished ok');

my $parse2 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $output_files[1],
    parse_outfile => $parse_files[1],

);

ok($parse2->execute, 'Parse2 finished ok, params option tested fine');

is(compare($parse_files[0], $parse_files[1]),0,'Parse1 and Parse2 are same');

my $parse3 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $output_files[2],
);

my $hsps = $parse3->execute;
isa_ok($hsps->[0],'Bio::Search::HSP::GenericHSP','Parse3 finished fine');
is($hsps->[0]->hit->seq_id, 'Contig18', 'Got right hit id');

my $parse4 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $test_dir.'/sample.blast.out',
    parse_outfile => $parse_files[2],
    percent_threshold => 99,
    length_threshold  => 1000,
);

my $hsps4 = $parse4->execute;
cmp_ok(scalar@$hsps4,'==',34, 'Parse4 finished ok, threshold filters are tested');
#isa_ok($hsps4->[0],'Bio::Search::HSP::GenericHSP','Got right Bio Perl HSP obj');



#$HeadURL$
#$Id$
