#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Temp;
use Test::More tests => 7;

BEGIN {
        use_ok('Genome::Model::Tools::WuBlast::Blastn');
        use_ok('Genome::Model::Tools::WuBlast::Parse');
}

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-WuBlast';
my $tmp_dir  = File::Temp::tempdir(
    "WuBlast_Blastn_XXXXXX", 
    TMPDIR => 1,
    CLEANUP => 1,
);

my @output_files = map{$tmp_dir."/blast$_.blast"}(1..3);

my $blast1 = Genome::Model::Tools::WuBlast::Blastn->create(
    database   => $test_dir.'/subject',
    query_file => $test_dir.'/query.fa', 
    output_file=> $output_files[0],
    E => '1e-150',
    M => 1, 
    N => -3,
    Q => 3,
    R => 3,
    W => 14,
    S2 => 2000,
);

ok($blast1->execute, 'Blast1 finished ok');

my $parse1 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $output_files[0],
);
my $hsps1 = $parse1->execute;

#my $out_file1 = $blast1->output_file;
#print "Blast1 out_put:$out_file1\n";

my $params_1 = 'Q=3 S2=2000 W=14 M=1 R=3 E=1e-150 N=-3';
my $params_2 = $params_1.' gapS2=10000 wordmask=seg lcmask -warnings -errors';

my $blast2 = Genome::Model::Tools::WuBlast::Blastn->create(
    database    => $test_dir.'/subject',
    query_file  => $test_dir.'/query.fa', 
    output_file => $output_files[1],
    params      => $params_1,
);

ok($blast2->execute, 'Blast2 finished ok, params option tested fine');

my $parse2 = Genome::Model::Tools::WuBlast::Parse->create(
    blast_outfile => $output_files[1],
);
my $hsps2 = $parse2->execute;

isa_ok($hsps1->[0], 'Bio::Search::HSP::GenericHSP', 'Got right HSP obj');
is_deeply($hsps1, $hsps2, 'Blast1 got the same results as Blast2');

my $blast3 = Genome::Model::Tools::WuBlast::Blastn->create(
    database    => $test_dir.'/subject',
    query_file  => $test_dir.'/query.fa', 
    params      => $params_2,
    output_file => $output_files[2],
);

ok($blast3->execute, 'Blast3 finished ok, other blast parameters tested fine');



#$HeadURL$
#$Id$
