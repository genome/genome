#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 4;

my $class = 'Genome::Notify';

use_ok($class);

class Genome::Notify::Notifier::Test {
    is => 'Genome::Notify::Notifier',
};

my $count;
my $last_notice;
sub Genome::Notify::Notifier::Test::notify {
    my $class = shift;
    my $notice = shift;

    $count++;
    $last_notice = $notice;

    isa_ok($notice, 'Genome::Notify::Notice');
    ok(!$notice->acknowledged, 'notice being notified is not acknowledged');
}

my $type = Genome::Notify::NoticeType->create(
    name => 'Test for Notify.t',
    default_notifier_class => 'Genome::Notify::Notifier::Test',
);

subtest 'setup a notice type' => sub {
    plan tests => 2;
    isa_ok($type, 'Genome::Notify::NoticeType');
    ok(!$type->__errors__, 'type is valid');
};

my %notice_info = (
    subject => UR::Value::Text->get('test'),
    type => $type,
    header => 'Test Notice!',
    body => 'This is a sentence in the test notice.',
    target => Genome::Sys->current_user,
);

subtest 'new notice' => sub {
    $count = 0;
    $class->notify(%notice_info);
    is($count, 1, 'notified once');
    $class->notify(%notice_info);
    is($count, 1, 'notice not repeated while not acknowledged');
    $last_notice->acknowledged(1);
    $class->notify(%notice_info);
    is($count, 2, 'recurrence of an acknowledge notice re-notifies');
    done_testing();
};

class Genome::Notify::Notifier::TestAlternate {
    is => 'Genome::Notify::Notifier',
};

my $alternate_count;
sub Genome::Notify::Notifier::TestAlternate::notify {
    my $class = shift;
    my $notice = shift;

    $alternate_count++;

    isa_ok($notice, 'Genome::Notify::Notice');
    ok(!$notice->acknowledged, 'notice being notified is not acknowledged');
}

subtest 'notice preference respected' => sub {
    my $preference = Genome::Notify::Preference->create(
        type => $type,
        notifier_class => 'Genome::Notify::Notifier::TestAlternate',
        target => Genome::Sys->current_user,
    );
    isa_ok($preference, 'Genome::Notify::Preference');
    ok(!$preference->__errors__, 'Preference has no errors');

    $count = 0;
    $alternate_count = 0;

    $notice_info{subject} = UR::Value::Text->get('preference test');

    $class->notify(%notice_info);
    is($count, 0, 'default notifier not triggered');
    is($alternate_count, 1, 'preferred notifier triggered instead');
    done_testing();
};
