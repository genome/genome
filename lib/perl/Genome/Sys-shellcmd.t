use strict;
use warnings;
use above 'Genome';
use Test::More;

test_redirect_stdout_stderr();
test_set_pipefail();

done_testing();

sub test_set_pipefail {
    my $cmd = 'echo $SHELLOPTS | grep pipefail';

    my $rv = Genome::Sys->shellcmd(cmd => $cmd);
    ok($rv, 'Default value of set_pipefail is True');

    $rv = Genome::Sys->shellcmd(cmd => $cmd, set_pipefail => 1);
    ok($rv, 'Can specify set_pipefail => 1');

    $rv = eval {Genome::Sys->shellcmd(cmd => $cmd, set_pipefail => 0)};
    ok(!$rv, 'Can specify set_pipefail => 0');
}

sub test_redirect_stdout_stderr {

    my $test_redirect = sub {
        my %params = @_;

        open my $savedout, '>&', \*STDOUT     or die "Can't dup STDOUT: $!";
        open my $savederr, '>&', \*STDERR     or die "Can't dup STDERR: $!";
        my $restore = UR::Util::on_destroy( sub {
            open(STDOUT, '>&', $savedout);
            open(STDERR, '>&', $savederr);
        });

        my($parentout, $childout, $parenterr, $childerr, $cmd);
        if ($params{parentout}) {
            $parentout = delete $params{parentout};
            #close(STDOUT);
            open(STDOUT, '>', $parentout);
            STDOUT->autoflush(1);
        }
        if ($params{parenterr}) {
            $parenterr = delete $params{parenterr};
            #close(STDERR);
            open(STDERR, '>', $parenterr);
            STDERR->autoflush(1);
        }

        return Genome::Sys->shellcmd(%params);
    };

    my $read_file = sub {
        my $name = shift;
        local($/);
        my $fh = IO::File->new($name, 'r');
        my $contents = <$fh>;
        return $contents;
    };

    my $do_dump_status_messages = Genome::Sys->dump_status_messages;
    Genome::Sys->dump_status_messages(0);
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 0;

    {
        my $parentout = File::Temp::tmpnam();
        my $childout = File::Temp::tmpnam();
        my $result = $test_redirect->(
                cmd => 'echo test',
                parentout => $parentout,
                redirect_stdout => $childout,
            );
        ok($result, 'Run echo with stdout redirected');
        ok(! (-s $parentout), 'Found no output in stdout of parent')
            or diag $read_file->($parentout);
        # two newlines, as the shellcmd will emit an extra newline
        is( $read_file->($childout), "test\n\n", 'redirected child stdout to a file');
        unlink($parentout, $childout);
    }

    {
        my $parenterr = File::Temp::tmpnam();
        my $childerr = File::Temp::tmpnam();
        my $result = $test_redirect->(
                cmd => 'echo test 1>&2',
                parenterr => $parenterr,
                redirect_stderr => $childerr,
            );
        ok($result, 'Run echo with stderr redirected');
        ok(! (-s $parenterr), 'Found no output in stderr of parent')
            or diag $read_file->($parenterr);
        is( $read_file->($childerr), "test\n", 'redirected child stderr to a file');
        unlink($parenterr, $childerr);
    }

    Genome::Sys->dump_status_messages($do_dump_status_messages);
}

