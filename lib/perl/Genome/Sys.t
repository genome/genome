#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Test::More;

use Genome::Sys;
use File::Temp;

sub mdir($) {
    system "mkdir -p $_[0]";
    ok(-d $_[0], "created directory $_[0]") or die "cannot continue!";
}

my $tmp = Genome::Sys->create_temp_directory("foo");
ok($tmp, "made temp directory $tmp");

my $tmp1 = $tmp . '/set1';
mdir($tmp1);
ok($tmp1, "made temp directory $tmp1");

my $tmp2 = $tmp . '/set2';
mdir($tmp2);
ok($tmp2, "made temp directory $tmp2");

$ENV{GENOME_DB} = join(":",$tmp1,$tmp2);

mdir($tmp1 . '/db1/1.0');
mdir($tmp1 . '/db1/2.1'); # the others are noise
mdir($tmp1 . '/db2/123');
mdir($tmp1 . '/db2/4');

my $ret = Genome::Sys->dbpath('db1','2.1');
is($ret, $tmp1 . '/db1/2.1', "path returns correctly");

mdir($tmp2 . '/db1/2.1'); # hidden by set 1

$ret = Genome::Sys->dbpath('db1','2.1');
is($ret, $tmp1 . '/db1/2.1', "path for db1 2.1 is the same as the last time because the new db is 2nd in the path");

rmdir $tmp1 . '/db1/2.1';
ok(! -d $tmp1 . '/db1/2.1', "removed the first database dir $tmp1/db1/2.1") or diag $!;

$ret = Genome::Sys->dbpath('db1','2.1');
is($ret, $tmp2 . '/db1/2.1', "path is the second db because the new db was removed") or diag $ret;

change_rollback_removes_symlink_for_create_symlink_and_log_change();

test_sudo_username();

test_file_operations();

done_testing();


sub change_rollback_removes_symlink_for_create_symlink_and_log_change {
    my $transaction = UR::Context::Transaction->begin();
    isa_ok($transaction, 'UR::Context::Transaction', 'transaction');

    my $object = UR::Value->get('foo');
    isa_ok($object, 'UR::Value', 'object');

    my $source = Genome::Sys->create_temp_directory();
    ok(-d $source, "source ($source) is a directory");

    my $destination_dir = Genome::Sys->create_temp_directory();
    ok(-d $destination_dir, "destination_dir ($destination_dir) is a directory");

    my $destination = $destination_dir . '/' . $object->id;

    Genome::Sys->create_symlink_and_log_change($object, $source, $destination);

    ok(-l $destination, "symlink created ($destination)");

    $transaction->rollback();

    ok(! -e $destination, "symlink destroyed in rollback");

    return 1;
}

sub test_sudo_username {
    no warnings qw(redefine);
    #Genome::Sys autoloaded here so it can be overridden
    my $username = Genome::Sys->username;

    {
        *Genome::Sys::cmd_output_who_dash_m = sub { return '' };
        local $ENV{SUDO_USER} = '';
        is(Genome::Sys->_sudo_username, '', 'sudo_username empty when not sudoed');
    }

    {
        *Genome::Sys::cmd_output_who_dash_m = sub { return '' };
        local $ENV{SUDO_USER} = "$username";
        is(Genome::Sys->_sudo_username, "$username", 'sudo_username detects based on SUDO_USER env var');
    }

    {
        *Genome::Sys::cmd_output_who_dash_m = sub { return "$username pt" };
        *Genome::Sys::username = sub { return "$username" };
        is(Genome::Sys->_sudo_username, '', 'sudo_username empty when not sudoed');
    }

    {
        *Genome::Sys::cmd_output_who_dash_m = sub { return "$username pt" };
        *Genome::Sys::username = sub { return 'not-user-name' };
        is(Genome::Sys->_sudo_username, "$username", 'sudo_username detects based on who -m');
    }

    use warnings qw(redefine);
}

sub test_file_operations {
    my $gzip_path = Genome::Sys->create_temp_file_path;
    ok ($gzip_path, "Got a tmp path");
    my $gzip_fh = Genome::Sys->open_gzip_file_for_writing($gzip_path);
    ok ($gzip_fh, "got a gzip fh");

    $gzip_fh->print("Testing");
    $gzip_fh->close;

    my $gzip_type = Genome::Sys->file_type($gzip_path);
    is($gzip_type, "gzip", "The file type is gzip");

    my $symlink_path = Genome::Sys->create_temp_file_path;
    ok($symlink_path, "Got first level symlink path");
    ok(Genome::Sys->create_symlink($gzip_path, $symlink_path), "Created first level symlink");

    my $second_symlink_path = Genome::Sys->create_temp_file_path;
    ok($second_symlink_path, "Got second level symlink path");
    ok(Genome::Sys->create_symlink($symlink_path, $second_symlink_path), "Created second level symlink");

    my $first_symlink_type = Genome::Sys->file_type($symlink_path);
    is($first_symlink_type, "gzip", "The symlink type is gzip");

    my $second_symlink_type = Genome::Sys->file_type($second_symlink_path);
    is($second_symlink_type, "gzip", "The second level symlink type is gzip");
}

