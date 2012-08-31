#!/usr/bin/env genome-perl
use strict;
use warnings;

use above "Genome";

use Test::More tests => 8;

use_ok("Genome::Sys::Email") or die "cannot contiue w/o the email module";

my $message = Genome::Sys::Email->get('apipe/2010-April/007022');
ok($message, 'got a mail object');

is($message->message_id,'007022', 'parsed out message id');
is($message->month,'2010-April', 'parsed out month of message');
is($message->list_name,'apipe', 'parsed out mailing list name');

ok($message->subject, 'loaded subject');
ok($message->body, 'loaded body');

my $non_message = Genome::Sys::Email->get('not an identifier');
ok(!$non_message, 'returned no object for invalid identifier');
