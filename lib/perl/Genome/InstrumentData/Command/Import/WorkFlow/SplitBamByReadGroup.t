#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require Genome::Utility::Test;
require File::Compare;
require File::Spec;
require File::Temp;
require List::MoreUtils;
require Sub::Install;
use Test::More;

my $class = 'Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup';
use_ok($class) or die;
my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', File::Spec->join('bam-rg-multi', 'v5')) or die;

my $autogenerate_new_object_id_uuid_sub = UR::Object::Type->can('autogenerate_new_object_id_uuid');
my $n = 0;
Sub::Install::reinstall_sub({ # to set the RG ID in the bam
        into => 'UR::Object::Type',
        as => 'autogenerate_new_object_id_uuid',
        code => sub{
            my @caller = caller;
            if ( $caller[0] eq $class ) {
                return ++$n x 32;
            }
            else {
                return $autogenerate_new_object_id_uuid_sub->();
            }
        },
    });

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $multi_rg_base_name = 'input.rg-multi.bam';
my $multi_rg_bam_path = File::Spec->join($tmp_dir, $multi_rg_base_name);
Genome::Sys->create_symlink(File::Spec->join($test_dir, $multi_rg_base_name), $multi_rg_bam_path);
ok(-s $multi_rg_bam_path, 'linked two read groups bam');
my $cmd = Genome::InstrumentData::Command::Import::WorkFlow::SplitBamByReadGroup->execute(bam_path => $multi_rg_bam_path);
ok($cmd->result, 'execute');

my @output_bam_paths = $cmd->output_bam_paths;
my @bam_basenames = (qw/ 2883581797.paired 2883581797.singleton 2883581798.paired 2883581798.singleton /);
is(@output_bam_paths, @bam_basenames, '4 read group bam paths');
for my $basename ( @bam_basenames ) {
    my $output_bam_path = File::Spec::->join($tmp_dir, 'input.rg-multi.'.$basename.'.bam');
    ok((List::MoreUtils::any { $_ eq $output_bam_path } @output_bam_paths), 'expected bam in output bams');
    ok(-s $output_bam_path, 'expectred bam path exists');
    my $expected_bam_path = File::Spec->join($test_dir, 'split-by-rg.'.$basename.'.bam');
    is(File::Compare::compare($output_bam_path, $expected_bam_path), 0, "expected $basename bam path matches");
}

#print "$tmp_dir\n"; <STDIN>;
done_testing();
