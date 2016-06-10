#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use File::Temp;
use Genome::Utility::Test 'compare_ok';
use Test::Exception;
use Test::More tests => 2;

my $pkg = 'Genome::Model::Tools::Sx::Trim::ByPosition';
use_ok($pkg) or die;

my $dir = File::Spec->join( Genome::Config::get('test_inputs'), 'Genome-Model-Tools-Sx', 'TrimByPosition');
my $tmp_dir = Genome::Sys->create_temp_directory;

subtest 'trim positions' => sub{
    plan tests => 4;

    my $positions_path = File::Spec->join($dir, 'trim-positions.dup.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Duplicate sequence id in trim positions\! seq2/, 'failed to load poistions with duplicates');

    my $positions_path = File::Spec->join($dir, 'trim-positions.invalid-stop.txt');
    throws_ok(sub{ $pkg->load_positions($positions_path); }, qr/Invalid stop position for seq3\! \?/, 'failed to load poistions with duplicates');

    $positions_path = File::Spec->join($dir, 'trim-positions.txt');
    my $positions = $pkg->load_positions($positions_path);
    ok($positions, 'load_positions');
    is_deeply(
        $positions,
        { seq1 => [[1, 'end']], seq2 => [[100, 'end']], seq3 => [[1, 100]], 'Contig1.1' => [[10, 20], [30,40]], },
        'trim_positions',
    );

};

done_testing();
