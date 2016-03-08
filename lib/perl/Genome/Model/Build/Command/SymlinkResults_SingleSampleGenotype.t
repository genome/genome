#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::SingleSampleGenotype;

use File::Spec;
use Test::More tests => 8;

my $class = 'Genome::Model::Build::Command::SymlinkResults';

use_ok($class, 'the class can be used');

my $build = setup_build();
my $destination = Genome::Sys->create_temp_directory();

my $cmd = $class->create(
    build => $build,
    destination => $destination,
);
isa_ok($cmd, $class, 'created a command');
ok($cmd->execute, 'command runs successfully');

my @destination_contents = glob(File::Spec->join($destination, '*'));
is(scalar(@destination_contents), 1, 'found a directory for the build');

my @build_contents = glob(File::Spec->join($destination_contents[0], '*'));
is(scalar(@build_contents), 3, 'created a directory for each type of results');
my $symlink_count = grep { -l $_ } @build_contents;
is($symlink_count, 2, 'found expected symlinks');
my ($dir) = grep { !-l $_ } @build_contents;
ok(-d $dir, 'found directory for HC results');

my @hc_contents = glob(File::Spec->join($dir, '*'));
my @intervals = test_intervals();
is(scalar(@hc_contents), scalar(@intervals), 'found expected number of HC results');


sub test_intervals {
    return (1..3, 'X', 'GL2000.1');
}

sub setup_build {
    my $test_model = Genome::Test::Factory::Model::SingleSampleGenotype->setup_object;
    my $test_build = Genome::Test::Factory::Build->setup_object(model_id => $test_model->id);

    my $id = -1;

    my $ss_result = Genome::InstrumentData::AlignmentResult::Merged::Speedseq->__define__(
        output_dir => Genome::Sys->create_temp_directory(),
        id => $id--,
    );
    Genome::SoftwareResult::User->create(
        software_result => $ss_result,
        user => $test_build,
        label => 'merged_alignment_result',
    );

    my $qc_result = Genome::Qc::Result->__define__(
        output_dir => Genome::Sys->create_temp_directory(),
        id => $id--,
    );
    Genome::SoftwareResult::User->create(
        software_result => $qc_result,
        user => $test_build,
        label => 'qc_result',
    );

    for my $interval (test_intervals()) {
        my $hc_result = Genome::Model::SingleSampleGenotype::Result::HaplotypeCaller->__define__(
            output_dir => Genome::Sys->create_temp_directory(),
            id => $id--,
        );
        $hc_result->add_input(name => 'intervals-0', value_id => $interval);
        Genome::SoftwareResult::User->create(
            software_result => $hc_result,
            user => $test_build,
            label => 'haplotype_caller_result',
        );
    }

    return $test_build;
}





