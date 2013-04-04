use strict;
use warnings;

use Test::Builder::Tester;
use Test::More tests => 10;
use UR;
use Sub::Install;

use_ok('Genome::Utility::Test', qw(command_execute_ok));

class My::Test::Command {
    is => 'Command::V2',
    has => [
        retval => { is => 'Integer' },
    ],
};
my $execute_cb = sub {};
Sub::Install::install_sub ({
    into => 'My::Test::Command',
    as => '_execute_body',
    code => sub {
        my $self = shift;
        $execute_cb->($self) if $execute_cb;
        return $self->retval;
    }
});

my $_command_execute_ok_parse_args = \&Genome::Utility::Test::_command_execute_ok_parse_args;

my($command, $message_config, $message);

eval {
    ($command, $message_config, $message)
        = $_command_execute_ok_parse_args->('My::Test::Command', {}, 'foo');
};
like($@,
    qr{Arg 1 to command_execute_ok\(\) must be an instance of a Command object},
    'Arg 1 must be a command object instance');

my $cmd = My::Test::Command->create( retval => 1 );
eval {
    ($command, $message_config, $message)
        = $_command_execute_ok_parse_args->($cmd, { error_messages => {} }, 'foo');
};
like($@,
    qr(Expected an arrayref for error_messages, but got HASH),
    'messages config must be an arrayref');

eval {
    ($command, $message_config, $message)
        = $_command_execute_ok_parse_args->($cmd, { error_messages => [ {} ] }, 'foo');
};
like($@,
    qr(Values for error_messages must be strings or Regexps),
    'Hashrefs not accepted as error_messages');


($command, $message_config, $message)
    = $_command_execute_ok_parse_args->($cmd);
is($command, $cmd, 'Parse with one arg returns command obj');
is_deeply($message_config, {}, 'Parse with one arg returns empty message_config');
is($message, 'execute My::Test::Command', 'Parse with one arg returns default test message');


($command, $message_config, $message)
    = $_command_execute_ok_parse_args->($cmd, 'hi there');
is($command, $cmd, 'Parse with two args returns command obj');
is_deeply($message_config, {}, 'Parse with two args returns empty message_config');
is($message, 'hi there', 'Parse with two args returns passed-in test message');


my $used_msg_config = {
    error_messages => [ 'one', qr(two) ],
    debug_messages => [ 'three' ],
    status_messages => undef,
};
($command, $message_config, $message)
    = $_command_execute_ok_parse_args->($cmd, $used_msg_config, 'ho there');
is($command, $cmd, 'Parse with three args returns command obj');
is_deeply($message_config, $used_msg_config, 'Parse with three args returns message config');
is($message, 'ho there', 'Parse with three args returns passed-in test message');


# test a successful execute and capture
$execute_cb = sub { shift->error_message('hi') };
$cmd = My::Test::Command->create( retval => 1 );
command_execute_ok($cmd,
        { error_messages => ['hi'] },
        'Match one error message');

# another success with multiple messages
$execute_cb = sub {
    my $self = shift;
    $self->error_message('hi'); $self->error_message('there');
    $self->status_message('foo');
};
$cmd = My::Test::Command->create( retval => 1 );
command_execute_ok($cmd,
        { error_messages => ['hi', 'there'], status_messages => ['foo'] },
        'Match multiple messages');

# try a regex
$execute_cb = sub { shift->error_message('hi there this is an error from the callback') };
$cmd = My::Test::Command->create( retval => 1 );
command_execute_ok($cmd,
        { error_messages => [qr(this is an error)] },
        'Match regex');


# Test a failed execute()
$execute_cb = sub {};
$cmd =  My::Test::Command->create( retval => 0 );
test_out('not ok 1 - blah (execute() returned false)');
test_fail(+1);
command_execute_ok($cmd, {}, 'blah');
test_test('execute() returning false is caught');


# test a message that does not match()
$execute_cb = sub { shift->error_message('hi') };
$cmd =  My::Test::Command->create( retval => 1 );
test_out('not ok 1 - blah');
test_fail(+2);
test_err(q(# For the 0th error_message, got 'hi' but expected 'bye'));
command_execute_ok($cmd, { error_messages => ['bye'] }, 'blah');
test_test('caught error_message that did not match');


# test expecting an extra message
$execute_cb = sub { shift->error_message('bye') };
$cmd =  My::Test::Command->create( retval => 1 );
test_out('not ok 1 - blah');
test_fail(+2);
test_err(q(# For the 1st error_message, got '' but expected 'another'));
command_execute_ok($cmd, { error_messages => ['bye', 'another'] }, 'blah');
test_test('expected error_message that did appear');


# test expecting no message but getting one
$execute_cb = sub { shift->error_message('bye') };
$cmd =  My::Test::Command->create( retval => 1 );
test_out('not ok 1 - blah');
test_fail(+2);
test_err(q(# For the 0th error_message, got 'bye' but expected ''));
command_execute_ok($cmd, { error_messages => [] }, 'blah');
test_test('got error message we did not expect');

