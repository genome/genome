#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Temp;
use Genome::Utility::Test 'compare_ok';
use Test::Exception;
use Test::More tests => 5;

my $pkg = 'Genome::Model::Tools::Sx::Trim::ByPosition';
use_ok($pkg) or die;

my $dir = File::Spec->join( Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Sx', 'TrimByPosition');
my $tmp_dir = Genome::Sys->create_temp_directory;

my %expected_positions = (
    TRIM_ALL => [],
    TRIM_LEFT => [ [1,60] ],
    TRIM_MIDDLE => [ [121,180] ],
    TRIM_RIGHT => [ [241,300] ],
    TRIM_LMR => [ [1,60],[121,180],[241,300] ],
    TRIM_MANY => [ [61,63],[121,123],[181,183] ],
    'TRIM_SPLIT.1' => [ [61,62] ],
    'TRIM_SPLIT.2' => [ [61,62] ],
    'Contig1.1' => [ [5,10] ],
    'Contig1.2' => [ [11,15] ],
);
subtest 'trim positions' => sub{
    plan tests => 4;

    my $positions_path = File::Spec->join($dir, 'trim-positions.dup.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Duplicate sequence id in trim positions\! seq2/, 'failed to load poistions with duplicates');

    $positions_path = File::Spec->join($dir, 'trim-positions.invalid-stop.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Invalid stop position for seq3\! \?/, 'failed to load poistions with duplicates');

    $positions_path = File::Spec->join($dir, 'trim-positions.txt');
    my $positions = $pkg->load_positions($positions_path);
    ok($positions, 'load_positions');
    is_deeply($positions, \%expected_positions, 'trim_positions');

};

subtest 'keep positions' => sub{
    plan tests => 7;

    my $trimmer = $pkg->create;
    ok($trimmer, 'create trimmer');

    my $positions = $trimmer->trim_positions(\%expected_positions);
    ok($trimmer->trim_positions, 'set positions');

    my $seq = { seq => 'A' x 300 };
    $seq->{id} = 'KEEP_ALL';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0, 300] ],
        'trim_positions for '.$seq->{id},
    );

    $seq->{id} = 'TRIM_ALL';
    is(
        $trimmer->keep_positions_for_sequence($seq),
        undef,
        'trim_positions for '.$seq->{id},
    );

    $seq->{id} = 'TRIM_MIDDLE';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,120], [180,120] ],
        'trim_positions for '.$seq->{id},
    );

    $seq->{id} = 'TRIM_MANY';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,60], [63,57], [123,57], [183,117], ],
        'trim_positions for '.$seq->{id},
    );

    $seq->{id} = 'Contig1.2';
    $seq->{orig_seq_id} = 'Scaffold1';
    is_deeply(
        $trimmer->keep_positions_for_sequence($seq),
        [ [0,10], [15,285] ],
        'trim_positions for '.$seq->{id},
    );

};

subtest 'trim sequence' => sub{
    plan tests => 8;

    my $trimmer = $pkg->create;
    ok($trimmer, 'create trimmer');

    my $positions = $trimmer->trim_positions(
        {
            TRIM_ALL => [],
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

subtest 'execute' => sub{
    plan tests => 3;

    my $positions_path = File::Spec->join($dir, 'trim-positions.txt');
    my $in_fasta = File::Spec->join($dir, 'input.fasta');
    my $out_fasta = File::Spec->join($tmp_dir, 'out.fasta');
    my $trimmer = $pkg->create(
        positions_path => $positions_path,
        input => $in_fasta,
        output => $out_fasta,
    );
    ok($trimmer, 'create trimmer');
    ok($trimmer->execute, 'execute');

    my $expected_fasta = File::Spec->join($dir, 'expected.fasta');
    compare_ok($out_fasta, $expected_fasta, 'out fasta matches');

};

done_testing();
