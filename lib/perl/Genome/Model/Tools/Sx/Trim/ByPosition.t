#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Temp;
use Genome::Utility::Test 'compare_ok';
use Test::Exception;
use Test::More tests => 4;

my $pkg = 'Genome::Model::Tools::Sx::Trim::ByPosition';
use_ok($pkg) or die;

my $dir = File::Spec->join( Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Sx', 'TrimByPosition');
my $tmp_dir = Genome::Sys->create_temp_directory;

subtest 'trim positions' => sub{
    plan tests => 4;

    my $positions_path = File::Spec->join($dir, 'trim-positions.dup.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Duplicate sequence id in trim positions\! seq2/, 'failed to load poistions with duplicates');

    $positions_path = File::Spec->join($dir, 'trim-positions.invalid-stop.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Invalid stop position for seq3\! \?/, 'failed to load poistions with duplicates');

    $positions_path = File::Spec->join($dir, 'trim-positions.txt');
    my $positions = $pkg->load_positions($positions_path);
    ok($positions, 'load_positions');
    is_deeply(
        $positions,
        { seq1 => 'ALL', seq2 => [[100, 'end']], seq3 => [[1, 100]], 'Contig1.1' => [[10, 20], [30,40]], },
        'trim_positions',
    );

};

subtest 'keep positions' => sub{
    plan tests => 7;

    my $trimmer = $pkg->create;
    ok($trimmer, 'create trimmer');

    my $positions = $trimmer->trim_positions(
        {
            TRIM_ALL => 'ALL',
            TRIM_ONE => [ [2,3] ],
            TRIM_SEVERAL => [ [2,4], [6,8], [11,13] ],
            'Contig1.1' => [ [3,6] ], # pcap
        },
    );
    ok($trimmer->trim_positions, 'set positions');

    my $seq = { seq => 'AAAAAAAAAAAAAAA' };
    $seq->{id} = 'KEEP_ALL';
    $seq->{orig_seq_id} = 'KEEP_ALL';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        'ALL',
        'trim_positions for KEEP_ALL',
    );

    $seq->{id} = 'TRIM_ALL';
    $seq->{orig_seq_id} = 'TRIM_ALL';
    is(
        $trimmer->keep_positions_for_sequence($seq),
        undef,
        'trim_positions for TRIM_ALL',
    );

    $seq->{id} = 'TRIM_ONE.1';
    $seq->{orig_seq_id} = 'TRIM_ONE';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,1], [3,12], ],
        'trim_positions for TRIM_ONE',
    );

    $seq->{id} = 'TRIM_SEVERAL.1';
    $seq->{orig_seq_id} = 'TRIM_SEVERAL';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,1], [4,1], [8,2], [13,2], ],
        'trim_positions for TRIM_SEVERAL',
    );

    $seq->{id} = 'Contig1.1';
    $seq->{orig_seq_id} = 'Scaffold1';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,2], [6,9], ],
        'trim_positions for Contig1.1',
    );

};

subtest 'trim sequence' => sub{
    plan tests => 8;

    my $trimmer = $pkg->create;
    ok($trimmer, 'create trimmer');

    my $positions = $trimmer->trim_positions(
        {
            TRIM_ALL => 'ALL',
            TRIM_SEVERAL => [ [2,4], [6,8], [11,13] ],
        },
    );
    ok($trimmer->trim_positions, 'set positions');

    my $seq = {
        seq => 'A' x 15,
        qual => '#' x 15,
    };
    $seq->{id} = 'KEEP_ALL';
    ok($trimmer->trim_sequence($seq), 'keep sequence');
    is_deeply(
        $seq,
        {
            id => 'KEEP_ALL',
            seq => 'A' x 15,
            qual => '#' x 15,
        },
        'sequence matches for KEEP_ALL',
    );

    $seq->{id} = 'TRIM_ALL';
    ok($trimmer->trim_sequence($seq), 'trim entire sequence');
    is_deeply(
        $seq,
        {
            id => 'TRIM_ALL',
            seq => '',
            qual => '',
        },
        'sequence matches for TRIM_ALL',
    );

    $seq->{id} = 'TRIM_SEVERAL';
    $seq->{seq} = 'A' x 15;
    delete $seq->{qual};
    #$seq->{qual} = '#' x 15;
    ok($trimmer->trim_sequence($seq), 'trim_sequence for TRIM_SEVERAL');
    is_deeply(
        $seq,
        {
            id => 'TRIM_SEVERAL',
            seq => 'A' x 6,
            #       qual => '#' x 6,
        },
        'sequence matches for TRIM_SEVERAL',
    );

};

done_testing();
