#!/usr/bin/env perl
use strict;
use warnings;
use above 'Genome';
use Genome::Utility::Test qw(command_execute_ok);

# This tests Genome::Command::WithSavedResults
# and also Genome::Command::BuildStepWrapper, which only
# wraps the former.

package Genome::TestCommand;

class Genome::TestCommand {
    is => 'Genome::Command::WithSavedResults',
    has_optional_param => [
        p1 => { is => 'Text' },
        p2 => { is => 'Number' },
        p3 => { is => 'Genome::Bar' },
        p4 => { is => 'Number', is_many => 1 },
        p5 => { is => 'Genome::Bar', is_many => 1 },
        p6 => { is => 'Boolean', default_value => 1 },
        result_version => { is => 'Integer', default_value => 1 },
    ],
    has_optional_input => [
        i1 => { is => 'Text' },
        i2 => { is => 'Number' },
        i3 => { is => 'Genome::Bar' },
        i4 => { is => 'Number', is_many => 1 },
        i5 => { is => 'Genome::Bar', is_many => 1 },
        output_dir => {},
    ],
    #has_metric => [
        # not yet supported
        #m1 => { is => 'Text' },
        #m2 => { is => 'Number' },
        #m3 => { is => 'Genome::Bar' },
        #m4 => { is => 'Number', is_many => 1 },
        #m5 => { is => 'Genome::Bar', is_many => 1 },
    #],
};

my $count = 0;
sub _execute_v1 {
    my $self = shift;
    $count++;
    $self->status_message("*** run $count: $self with p1 " . $self->p1);
    return 1;
}

package main;
use Test::More tests => 13;

my $result_meta = Genome::TestCommand::Result->__meta__;
ok($result_meta, "got result meta for new class");

my $wrap_meta = Genome::TestCommand::BuildStepWrapper->__meta__;
ok($wrap_meta, "successfully generated a wrapper for generating software results within builds");

my $build;
{
    local $SIG{__WARN__} = sub {
        my $msg = shift;
        die "Unexpected warning message: $msg" unless($msg =~ m/^exists:/);
    };
    $build = Genome::Model::Build->get(135315483);
}
ok($build, "got test build");
print join("\n",Genome::Model::Build->error_messages),"\n";

# temp change data dir to test location
my $dir1 = Genome::Sys->create_temp_directory();
$build->data_directory($dir1);

my $meta1 = Genome::TestCommand::BuildStepWrapper->__meta__;
ok($meta1, "got BuildStepWrapper for command class");

my $dir2 = Genome::Sys->create_temp_directory();

Genome::TestCommand->dump_status_messages(0);
Genome::TestCommand::BuildStepWrapper->dump_status_messages(0);
Genome::SoftwareResult::Stageable->dump_status_messages(0);
my $wrapper_result1 = Genome::TestCommand::BuildStepWrapper->execute(
    p1 => "P1", 
    i1 => "I1", 
    output_dir => $dir2,
    wrapper_build => $build, 
    wrapper_build_label => "test_label", 
    result_version => 1,
);
is($count, 1, "the underlying command has run one time");

Genome::TestCommand::Result->dump_error_messages(0);
Genome::TestCommand::Result->queue_error_messages(1);
my $sr1 = Genome::TestCommand::Result->get(p1 => "P1", i1 => "I1", result_version => 1);
ok($sr1, "got result for wrapper");
my @error_messages = Genome::TestCommand::Result->error_messages();
is(scalar(@error_messages), 1, 'Expected 1 error messages');
like($error_messages[0],
    qr/^Calling get on SoftwareResult \(unless getting by id\) is slow and possibly incorrect\./,
    'Expected error text');

my @u1 = $sr1->users();
is(scalar(@u1), 1, "got one 'user' of the result");

is($u1[0]->label, "test_label", "the label is correct linking the build to the SR");
is($u1[0]->user, $build, "user of the result on that link is the build");
is($u1[0]->software_result, $sr1, "the SR on the link is the one we created");

my $dir3 = Genome::Sys->create_temp_directory();
$build->data_directory($dir3);

my $dir4 = Genome::Sys->create_temp_directory();

my $wrapper_result2 = Genome::TestCommand::BuildStepWrapper->execute(
    p1 => "P1", 
    i1 => "I1", 
    output_dir => $dir4,
    wrapper_build => $build, 
    wrapper_build_label => "test_label", 
    result_version => 1,
);
is($count, 1, "the underlying command has still run just one time");


1;

