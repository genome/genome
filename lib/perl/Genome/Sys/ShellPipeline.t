#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use Test::Exception;
use File::Slurp qw(read_file);
use above 'Genome';

my $pkg = 'Genome::Sys::ShellPipeline';
use_ok($pkg);

my $tmpdir = Genome::Sys->create_temp_directory;

sub make_subtest {
    my ($name, @args) = @_;

    $tmpdir = Genome::Sys->create_temp_directory;

    subtest $name => @args;
}

make_subtest "success no capture" => sub {
    my $obj = $pkg->create(
        pipe_commands => [
            qq{$^X -e 'print "hello\\n"'},
            qq{$^X -e 'my \$line = <>; print STDERR uc \$line; print \$line;'}
            ],
        );

    ok($obj, "created command");
    ok($obj->execute, "executed pipeline");
};

make_subtest "success capture" => sub {
    my $tmp_stdout = File::Spec->catfile($tmpdir, "stdout.txt");
    my $tmp_stderr = File::Spec->catfile($tmpdir, "stderr.txt");

    my $obj = $pkg->create(
        pipe_commands => [
            qq{$^X -e 'print "hello\\n"'},
            qq{$^X -e 'my \$line = <>; print STDERR uc \$line; print \$line;'}
            ],
        redirects => "> $tmp_stdout 2> $tmp_stderr",
        );

    ok($obj, "created command");
    ok($obj->execute, "executed pipeline");

    is(read_file($tmp_stdout), "hello\n", "stdout captured");
    is(read_file($tmp_stderr), "HELLO\n", "stderr captured");
};

make_subtest "success capture w/pre post" => sub {
    my $tmp_pre = File::Spec->catfile($tmpdir, "pre.txt");
    my $tmp_post = File::Spec->catfile($tmpdir, "post.txt");
    my $tmp_stdout = File::Spec->catfile($tmpdir, "stdout.txt");
    my $tmp_stderr = File::Spec->catfile($tmpdir, "stderr.txt");

    my $obj = $pkg->create(
        pre_commands => [
            qq{echo "no pre capture"},
            qq{touch $tmp_pre},
        ],
        post_commands => [
            qq{echo "no post capture"},
            qq{touch $tmp_post},
        ],

        pipe_commands => [
            qq{$^X -e 'print "hello\\n"'},
            qq{$^X -e 'my \$line = <>; print STDERR uc \$line; print \$line;'}
            ],

        redirects => "> $tmp_stdout 2> $tmp_stderr",
        );

    ok($obj, "created command");
    ok($obj->execute, "executed pipeline");

    is(read_file($tmp_stdout), "hello\n", "stdout captured");
    is(read_file($tmp_stderr), "HELLO\n", "stderr captured");

    ok(-f $tmp_pre, "pre commands were executed");
    ok(-f $tmp_post, "post commands were executed");
};

make_subtest "pre failure ignored by default" => sub {
    my $tmp_post = File::Spec->catfile($tmpdir, "post.txt");
    my $tmp_stdout = File::Spec->catfile($tmpdir, "stdout.txt");
    my $tmp_stderr = File::Spec->catfile($tmpdir, "stderr.txt");

    my $obj = $pkg->create(
        pre_commands => [make_return_value_command(1)],
        post_commands => [
            qq{echo "no post capture"},
            qq{touch $tmp_post},
        ],

        pipe_commands => [
            qq{$^X -e 'print "hello\\n"'},
            qq{$^X -e 'my \$line = <>; print STDERR uc \$line; print \$line;'}
            ],

        redirects => "> $tmp_stdout 2> $tmp_stderr",
        );

    ok($obj, "created command");
    ok($obj->execute, "executed pipeline");

    is(read_file($tmp_stdout), "hello\n", "stdout captured");
    is(read_file($tmp_stderr), "HELLO\n", "stderr captured");

    ok(-f $tmp_post, "post commands were executed");
};

make_subtest "pre failure can be caught with bash" => sub {
    my $tmp_post = File::Spec->catfile($tmpdir, "post.txt");
    my $tmp_stdout = File::Spec->catfile($tmpdir, "stdout.txt");
    my $tmp_stderr = File::Spec->catfile($tmpdir, "stderr.txt");

    my $obj = $pkg->create(
        pre_commands => ["set -e", make_return_value_command(1)],
        post_commands => [
            qq{echo "no post capture"},
            qq{touch $tmp_post},
        ],

        pipe_commands => [
            qq{$^X -e 'print "hello\\n"'},
            qq{$^X -e 'my \$line = <>; print STDERR uc \$line; print \$line;'}
            ],

        redirects => "> $tmp_stdout 2> $tmp_stderr",
        );

    ok($obj, "created command");
    dies_ok {$obj->execute} "command failed as expected";
};

make_subtest "failure detection" => sub {
    my $tmp_pre = File::Spec->catfile($tmpdir, "pre.txt");
    my $tmp_post = File::Spec->catfile($tmpdir, "post.txt");

    my $ok_cmd = qq{$^X -e 'exit 0;'};
    my $fail_cmd = qq{$^X -e 'exit 1;'};

    # Test all combinations of 0/1 return codes $n_commands processes.
    my $n_commands = 4; # 2^4 = 16 tests total
    for my $i (0..(2**$n_commands - 1)) {
        my @rvs = split //, sprintf("%0${n_commands}b", $i);
        my @cmds = map {make_return_value_command($_)} @rvs;

        my $cmd = $pkg->create(
            pre_commands => [qq{touch $tmp_pre}],
            post_commands => [qq{touch $tmp_post}],
            pipe_commands => \@cmds
            );

        ok($cmd, "created command");

        unlink $tmp_post;

        if ($i != 0) {
            throws_ok {$cmd->execute} qr/The following commands crashed:/, "command crashed as expected";
            ok(!-e $tmp_post, "post command was not executed");
        }
        else {
            ok($cmd->execute, "executed command");
            ok(!$@, "command succeeded as expected");
            ok(-e $tmp_post, "post command was executed");
        }

        ok(-f $tmp_pre, "pre commands were executed");
        is_deeply($cmd->return_codes, \@rvs, "return codes match");
    }
};

done_testing();

sub make_return_value_command {
    my $rv = shift;
    return qq{$^X -e "exit $rv;"};
}
