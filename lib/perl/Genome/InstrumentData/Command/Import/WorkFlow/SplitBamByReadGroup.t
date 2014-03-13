#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Temp;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup') or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import') or die;

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $multi_rg_base_name = 'input.rg-multi';
my $multi_rg_bam_path = $tmp_dir.'/'.$multi_rg_base_name.'.bam';
Genome::Sys->create_symlink($test_dir.'/bam-rg-multi/v3/'.$multi_rg_base_name.'.bam', $multi_rg_bam_path);
ok(-s $multi_rg_bam_path, 'linked two read groups bam');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $multi_rg_bam_path);
ok($cmd, 'execute');
my @read_group_bam_paths = $cmd->read_group_bam_paths;
is(@read_group_bam_paths, 4, '4 read group bam paths');
is_deeply(
    [sort @read_group_bam_paths], 
    [sort map { $tmp_dir.'/'.$multi_rg_base_name.'.'.$_.'.bam' } (qw/ 2883581797.paired 2883581797.singleton 2883581798.paired 2883581798.singleton /)],
    '4 read groups bams paths');
is(grep({ -s } @read_group_bam_paths), 4, 'read group bam paths exist');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
