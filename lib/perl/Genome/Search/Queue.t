#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;

use above "Genome";
use_ok('Genome::Search::Queue') || die;

require Genome::Search;
my $orig_is_indexable = \&Genome::Search::is_indexable;
my $text_is_indexable = sub {
    my ($class, $object) = @_;
    if ($object->isa('UR::Value::Text')) {
        return 1;
    }
    else {
        return;
    }
};

test_create_missing_subject();
test_create_missing_timestamp();
test_create_existing_subject();
test_create_non_indexable_subject();
test_priority_sorting();

done_testing();

sub test_create_missing_subject {
    my $tx = UR::Context::Transaction->begin();

    my $index_queue = eval {
        Genome::Search::Queue->create();
    };
    my $error = $@;

    is($index_queue, undef, 'failed to create index_queue when missing subject');
    like($error, qr/subject/, 'error mentions subject');

    $tx->rollback();
}

sub test_create_missing_timestamp {
    my $tx = UR::Context::Transaction->begin();

    *Genome::Search::is_indexable = $text_is_indexable;

    my $subject = UR::Value::Text->get('Hello, world.');
    my $index_queue = Genome::Search::Queue->create(
        subject_id => $subject->id,
        subject_class => $subject->class,
    );

    isa_ok($index_queue, 'UR::Object', 'create returned an object');
    ok($index_queue->timestamp, 'timestamp was added');

    *Genome::Search::is_indexable = $orig_is_indexable;

    $tx->rollback();
}

sub test_create_existing_subject {
    my $tx = UR::Context::Transaction->begin();

    *Genome::Search::is_indexable = $text_is_indexable;

    my $subject = UR::Value::Text->get('Hello, world.');
    my $index_queue = Genome::Search::Queue->create(
        subject_id => $subject->id,
        subject_class => $subject->class,
    );

    isa_ok($index_queue, 'UR::Object', 'create returned an object');

    my $index_queue_2 = Genome::Search::Queue->create(
        subject_id => $subject->id,
        subject_class => $subject->class,
    );
    isa_ok($index_queue_2, 'UR::Object', 'create returned an object');

    isnt($index_queue_2, $index_queue, 'new index_queue_2 is different than index_queue');

    *Genome::Search::is_indexable = $orig_is_indexable;

    $tx->rollback();
}

sub test_create_non_indexable_subject {
    my $tx = UR::Context::Transaction->begin();

    my $subject = UR::Value::Text->get('Hello, world.');
    my $index_queue = eval {
        Genome::Search::Queue->create(
            subject_id => $subject->id,
            subject_class => $subject->class,
        );
    };
    my $error = $@;

    is($index_queue, undef, 'failed to create index_queue when subject is not indexable');
    like($error, qr/indexable/, 'error mentions indexable');

    $tx->rollback();
}

sub test_priority_sorting {
    my $tx = UR::Context::Transaction->begin();

    *Genome::Search::is_indexable = $text_is_indexable;

    # purposely out of order so that timestamps won't sort in the same order as priority
    my @subject_ids = ('Thing 9', 'Thing', 'Thing 5', 'Thing 1', 'Thing 0', 'Thing 2');

    # undef = 0, numerically ascending
    my @sorted_subject_ids = ('Thing 0', 'Thing 1', 'Thing 2', 'Thing 5', 'Thing 9', 'Thing');

    for my $subject_id (@subject_ids) {
        my ($priority) = $subject_id =~ /(\d)$/;
        my $subject = UR::Value::Text->get($subject_id);
        my $index_queue = Genome::Search::Queue->create(
            subject_id => $subject->id,
            subject_class => $subject->class,
            priority => $priority,
        );
        isa_ok($index_queue, 'Genome::Search::Queue', "created queue object for $subject_id");
        sleep 1;
    }

    my @priority_sorted_queue = Genome::Search::Queue->get(subject_id => \@subject_ids, -order_by => 'priority');
    ok(@priority_sorted_queue, 'got priority_sorted_queue') || return;
    my @priority_sorted_queue_subject_ids = map { $_->subject_id } @priority_sorted_queue;
    is_deeply(\@priority_sorted_queue_subject_ids, \@sorted_subject_ids, 'ordering by priority matches expected results')
        or diag(Data::Dumper::Dumper(\@priority_sorted_queue_subject_ids, \@sorted_subject_ids));

    my @timestamp_sorted_queue = Genome::Search::Queue->get(subject_id => \@subject_ids, -order_by => 'timestamp');
    ok(@priority_sorted_queue, 'got timestamp_sorted_queue') || return;
    my @timestamp_sorted_queue_subject_ids = map { $_->subject_id } @timestamp_sorted_queue;
    isnt(join('', @timestamp_sorted_queue_subject_ids), join('', @sorted_subject_ids), 'ordering by timestamp does not match priorty sort');
    isnt(join('', @timestamp_sorted_queue_subject_ids), join('', @priority_sorted_queue_subject_ids), 'ordering by timestamp does not match priorty sort');

    $tx->rollback();
}
