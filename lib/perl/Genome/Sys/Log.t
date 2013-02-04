#!/usr/bin/env perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 294;
use JSON;

#
# NOTE: this tests only the experimental/optional JSON logging 
# which occurs when GENOME_SYS_LOG_DETAIL is set
#

# this will cause us to cut out just before sending into the syslogger
# but will test the rest of our decision making and processing
$Genome::Sys::Log::test_syslog = 1;

# this turns on JSON logging
$ENV{GENOME_SYS_LOG_DETAIL} = 1;

# the code which logs is in F2::f2, called by F1::f1, called below in namespace F0
require __FILE__ . ".d/F1.pm";
require __FILE__ . ".d/F2.pm";

# try all possible permutations of message emission and specified log level, 
# including having logging turned off 
my $n = 0;
my $type_n = 0;
for my $type (qw/debug status warning error/) {
    $type_n++;

    my $level = $type;
    if ($level eq 'status') {
        $level = 'info';
    }
    elsif ($level eq 'warning') {
        $level = 'warn';
    }

    my $level_n = 0;
    for my $min_level (qw/debug info warn error/,'') {
        $level_n++;

        $ENV{GENOME_SYS_LOG_LEVEL} = $min_level;
       
        for my $form (qw/text json/) {
            $n++;
            
            $Genome::Sys::Log::last_syslog_retval = undef;
            $Genome::Sys::Log::last_syslog_level = undef;
            $Genome::Sys::Log::last_syslog_message = undef;

            my $msg;
            if ($form eq 'text') {
                $msg = "test text $n $type";
            }
            else {
                $msg =  encode_json({ foo => "test json $n $type", n => $n, tt => $type });
            }

            my $call_line = __LINE__ + 2;
            package F0;
            my $retval = F1::f1($type, $msg);
            package main;

            if ($type_n >= $level_n) {
                is($Genome::Sys::Log::last_syslog_retval, 1, "message $n for type $type got retval $retval");
                is($Genome::Sys::Log::last_syslog_level, $level, "message $n for type $type is level $level");

                my $data = ($Genome::Sys::Log::last_syslog_message ? decode_json($Genome::Sys::Log::last_syslog_message) : {});
                is($data->{type}, $level, "message $n has level $level as expected");
                is($data->{u}, $ENV{USER}, "message $n has user set to $ENV{USER} as expected");
                is($data->{p}, $$, "message $n has pid set to $$ as expected"); 
                ok($data->{f} =~ /\/F2.pm$/, "message $n has correct filename $data->{f}");
                is($data->{l}, 7, "message $n has expected line 7");
                is($data->{m}, "F2::f2", "message $n has expected method F2::f2");

                my $e1 = "$call_line|F1::f1,5|F2::f2,7";

                is(index($data->{c}, $e1)+length($e1),length($data->{c}), "call stack string has expected value '$e1' (possibly with some prefix)");
               
                if ($form eq 'json') {
                    is($data->{n}, $n, "message $n has test field n with the correct value");
                    is($data->{tt}, $type, "message $n has tes field tt with the correct value");
                    is($data->{foo}, "test json $n $type");
                }
                else {
                    is($data->{msg}, $msg, "message $n for type $type is $msg");
                }
            }
            else {
                is($Genome::Sys::Log::last_syslog_retval, undef, "no logging -> no retval");
                is($Genome::Sys::Log::last_syslog_level, undef, "no loggging -> no log level");
                is($Genome::Sys::Log::last_syslog_message, undef, "no logging -> no message");
            }
        }
    }
}

#
# test the special override so that Genome::Sys messages are linked to their caller not the Genome::Sys line
#

$ENV{GENOME_SYS_LOG_LEVEL} = 'debug';
my $expected_line1;
my $expected_line2;
sub F0::f0 {
    $expected_line2 = __LINE__ + 1;
    Genome::Sys->shellcmd(cmd => 'sleep 1');
}
$expected_line1 = __LINE__ + 1;
F0::f0();

is($Genome::Sys::Log::last_syslog_retval, 1, "log retval is true");
is($Genome::Sys::Log::last_syslog_level, 'debug', "the last message is of type debug");

my $data = ($Genome::Sys::Log::last_syslog_message ? decode_json($Genome::Sys::Log::last_syslog_message) : {});
is($data->{type}, 'debug', "");
is($data->{u}, $ENV{USER}, "");
is($data->{p}, $$, ""); 
ok($data->{f} =~ /Log.t/, "file is this file, not the Genome::Sys file: got $data->{f}");
is($data->{l}, $expected_line2, "the line number is the line of our sys call, not the line within which emits the message");
is($data->{m}, "F0::f0", "method is F0::f0");
ok(index($data->{c},"$expected_line1|F0::f0,$expected_line2") == 0, "call stack string has expected start value '$expected_line1|F0::f0,$expected_line2'");
ok(defined($data->{t1}), "t1 has a value");
ok(defined($data->{t2}), "t2 has a value");
ok(defined($data->{elapsed}), "elapsed has a value");
ok($data->{t2} >= $data->{t1}, "t2 is not before t1");
ok($data->{t2} - $data->{t1} == $data->{elapsed}, "elapsed time is correct");

