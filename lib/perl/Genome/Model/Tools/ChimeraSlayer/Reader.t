#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require IO::String;
use Test::More;

use_ok('Genome::Model::Tools::ChimeraSlayer::Reader') or die;

my $version = 1;
my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ChimeraSlayer/v'.$version;
my $cpc_file = $test_data_dir.'/chims1.NAST.CPS.CPC';

# 0  ChimeraSlayer
# 1  chimera_X92624|S000015368_d10.86_AF282889|S000390572_d11.18_nc4580_ec939
# 2  7000004131495956
# 3  S000469847
# 4  1.0360
# 5  99.35
# 6  100
# 7  0.9354
# 8  89.70
# 9  0
# 10  YES
# 11 NAST:4595-4596
# 12 ECO:941-942

my $reader = Genome::Model::Tools::ChimeraSlayer::Reader->create(
    input => $cpc_file,
);
ok($reader, 'create');
my @expected_chimeras = _expected_chimeras();
my $cnt = 0;
while ( my $chimera = $reader->read ) {
    ok($chimera, "got chimera $cnt");
    #print Dumper($chimera);
    is_deeply(
        $chimera,
        $expected_chimeras[$cnt],
        "chimera $cnt matches",
    );
    $cnt++;
}

# fails
my $io = IO::String->new();
$io->print("ChimeraSlayer\tblah\n");
$io->print("ChimeraSlayer\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n");
$io->seek(0, 0);
$reader = Genome::Model::Tools::ChimeraSlayer::Reader->create(
    input => $io,
);
my $rv = eval{ $reader->read; };
ok(!$rv, 'failed to read line does not have enough fields');
like($@, qr/Malformed chimera slayer line! Got 1 fields, but expected 10 or 12 from chimera slayer line:/, 'error message matches');
$rv = eval{ $reader->read; };
ok(!$rv, 'failed to read line where verdict was not YES/NO');
like($@, qr/Malformed chimera slayer line! Verdict is expected to be YES, NO or UNKNOWN but is 10./, 'error message matches');

done_testing();


###

sub _expected_chimeras {
    return (
        {
            'percent_identity_right_A_left_B' => '89.70',
            'percent_identity_left_A_right_B' => '99.35',
            'parent_A' => '7000004131495956',
            'verdict' => 'YES',
            'parent_B' => 'S000469847',
            'est_breakpoint_ECOLI' => 'ECO:941-942',
            'est_breakpoint_NAST' => 'NAST:4595-4596',
            'divergence_ratio_right_A_left_B' => '0.9354',
            'id' => 'chimera_X92624|S000015368_d10.86_AF282889|S000390572_d11.18_nc4580_ec939',
            'divergence_ratio_left_A_right_B' => '1.0360',
            'confidence_left_A_right_B' => '100',
            'confidence_right_A_left_B' => '0'
        },
        {
            'percent_identity_right_A_left_B' => '89.26',
            'percent_identity_left_A_right_B' => '99.37',
            'parent_A' => 'S000414716',
            'verdict' => 'YES',
            'parent_B' => 'S000276014',
            'est_breakpoint_ECOLI' => 'ECO:918-919',
            'est_breakpoint_NAST' => 'NAST:4552-4553',
            'divergence_ratio_right_A_left_B' => '0.9343',
            'id' => 'chimera_L37602|S000414716_d11.44_AJ535638|S000276014_d11.46_nc4544_ec918',
            'divergence_ratio_left_A_right_B' => '1.0401',
            'confidence_left_A_right_B' => '100',
            'confidence_right_A_left_B' => '0'
        },
        {
            'percent_identity_right_A_left_B' => '99.65',
            'percent_identity_left_A_right_B' => '89.20',
            'parent_A' => 'S000270785',
            'verdict' => 'YES',
            'parent_B' => 'S000020199',
            'est_breakpoint_ECOLI' => 'ECO:726-727',
            'est_breakpoint_NAST' => 'NAST:3876-3881',
            'divergence_ratio_right_A_left_B' => '1.0538',
            'id' => 'chimera_X87340|S000020199_d10.88_AB184094|S000651818_d10.81_nc3822_ec703',
            'divergence_ratio_left_A_right_B' => '0.9433',
            'confidence_left_A_right_B' => '0',
            'confidence_right_A_left_B' => '100'
        },
    );
}
