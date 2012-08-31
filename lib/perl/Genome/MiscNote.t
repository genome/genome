#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

sudo_username_is_detected();
sudo_username_isnt_detected();
no_auto_truncate_body_text_does_not_truncates_long_text();
auto_truncate_body_text_truncates_long_text();
auto_truncate_body_text_truncates_long_text_with_sudo_username();

done_testing();

sub miscnote_body_text {
    my $sudo_username = shift;
    my $body_text = shift;

    no warnings qw(redefine);
    *Genome::Sys::sudo_username = sub { return $sudo_username };
    use warnings qw(redefine);

    my $subject = UR::Value->get('subject');
    isa_ok($subject, 'UR::Value', 'subject');

    my $note = Genome::MiscNote->create(
        subject => $subject,
        header_text => 'Test Note',
        body_text => $body_text,
    );
    isa_ok($note, 'Genome::MiscNote', 'note');

    return $note->body_text;
}

sub sudo_username_isnt_detected {
    my $message = 'Sample note message.';
    my $body_text = miscnote_body_text('', $message);
    is($body_text, $message, 'note didnt prepend sudo_username');
}

sub sudo_username_is_detected {
    my $message = 'Sample note message.';
    my $sudo_username = 'sample-sudo-username';
    my $body_text = miscnote_body_text($sudo_username, $message);
    is($body_text, $sudo_username . ' is running as ' . Genome::Sys->username . ". $message", 'note prepended with sudo_username');
}

sub no_auto_truncate_body_text_does_not_truncates_long_text {
    my $max_length = Genome::MiscNote->__meta__->property('body_text')->data_length;
    my $body_text = "A"x($max_length + 1);

    my $subject = UR::Value->get('subject');
    isa_ok($subject, 'UR::Value', 'subject');

    my $note = Genome::MiscNote->create(
        subject => $subject,
        header_text => 'Test Note',
        body_text => $body_text,
        auto_truncate_body_text => 0,
    );
    isa_ok($note, 'Genome::MiscNote', 'note');
    is(length($note->body_text), ($max_length + 1), "body_text was not truncated to $max_length");
}
sub auto_truncate_body_text_truncates_long_text {
    my $max_length = Genome::MiscNote->__meta__->property('body_text')->data_length;
    my $body_text = "A"x($max_length + 1);

    my $subject = UR::Value->get('subject');
    isa_ok($subject, 'UR::Value', 'subject');

    my $note = Genome::MiscNote->create(
        subject => $subject,
        header_text => 'Test Note',
        body_text => $body_text,
        auto_truncate_body_text => 1,
    );
    isa_ok($note, 'Genome::MiscNote', 'note');
    is(length($note->body_text), $max_length, "body_text was truncated to $max_length");
}
sub auto_truncate_body_text_truncates_long_text_with_sudo_username {
    my $max_length = Genome::MiscNote->__meta__->property('body_text')->data_length;
    my $body_text = "A"x($max_length + 1);
    my $sudo_username = 'sample-sudo-username';
    my $subject = UR::Value->get('subject');
    isa_ok($subject, 'UR::Value', 'subject');

    no warnings qw(redefine);
    *Genome::Sys::sudo_username = sub { return $sudo_username };
    use warnings qw(redefine);

    my $note = Genome::MiscNote->create(
        subject => $subject,
        header_text => 'Test Note',
        body_text => $body_text,
        auto_truncate_body_text => 1,
    );
    isa_ok($note, 'Genome::MiscNote', 'note');
    is(length($note->body_text), $max_length, "body_text was truncated to $max_length");
}
