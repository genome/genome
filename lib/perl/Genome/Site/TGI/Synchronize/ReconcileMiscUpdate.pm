package Genome::Site::TGI::Synchronize::ReconcileMiscUpdate;

use strict;
use warnings;

use Genome;

use Date::Format;
use Date::Parse;

class Genome::Site::TGI::Synchronize::ReconcileMiscUpdate {
    is => 'Command::V2',
    has => [
        start_from => { 
            is => 'Text', 
            doc => 'Date to start from to get updates to reconcile. Format is YYYY-MM-DD HH:MM::SS.',
        },
    ],
    has_optional_transient => [
        _misc_updates => { is => 'Array', },
        _updates => { is => 'Array', },
        _indels => { is => 'Array', },
        _stop_at => { is => 'Text', },
        stats => { is => 'Hash', default_value => {}, },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    $self->_stop_at( Date::Format::time2str("%Y-%m-%d %X", time()) ) if not $self->_stop_at; # set to now, unless overridden

    my $date = $self->start_from;
    my $time = Date::Parse::str2time($date);
    if ( not $time ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ start_from /],
            desc => 'Invalid date format => '.$date,
        );
    }

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Reconcile Misc Updates...');

    my $load_ok = $self->_load_misc_updates;
    return if not $load_ok;

    my $execute_updates = $self->_execute_updates;
    return if not $execute_updates;

    my $execute_indels = $self->_execute_indels;
    return if not $execute_indels;

    #$self->_execute_report;

    $self->status_message('Reconcile Misc Updates...done');
    return 1;
}

sub _load_misc_updates {
    my $self = shift;

    return $self->_misc_updates if $self->_misc_updates;

    $self->status_message('Load misc updates...');
    $self->status_message('Start date: '.$self->start_from);
    $self->status_message('Stop date:  '.$self->_stop_at);
    my @misc_updates = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->get(
        'edit_date between' => [ $self->start_from, $self->_stop_at ],
        is_reconciled => 0,
        '-order_by' => 'edit_date',
    );
    $self->status_message('Total found: '.@misc_updates);

    $self->status_message('Load misc updates...done');
    return $self->_misc_updates(\@misc_updates);
}

sub _execute_updates {
    my $self = shift;
    $self->status_message('Execute UPDATES...');

    my @misc_updates = grep { $_->description eq 'UPDATE' } @{$self->_misc_updates};
    $self->status_message('UPDATES found: '.@misc_updates);
    for my $misc_update ( @misc_updates ) {
        $misc_update->perform_update;
    }
    $self->_updates(\@misc_updates);

    $self->status_message('Execute UPDATES...done');
    return 1;
}

sub _execute_indels {
    my $self = shift;
    $self->status_message('Execute INDELS...');

    my @misc_updates = grep { $_->description eq 'INSERT' or $_->description eq 'DELETE' }  @{$self->_misc_updates};
    my (%subject_attr_misc_updates, @order);
    foreach my $misc_update ( @misc_updates ) {
        my $subject_attr_misc_update = Genome::Site::TGI::Synchronize::Classes::MiscUpdate::SubjectAttribute->get_or_create_from_misc_updates($misc_update);
        if ( not $subject_attr_misc_updates{ $subject_attr_misc_update->id } ) {
            $subject_attr_misc_updates{ $subject_attr_misc_update->id } = $subject_attr_misc_update;
            push @order, $subject_attr_misc_update->id;
        }
    }
    $self->status_message('INDELS found: '.@order);

    my @misc_update_indels;
    for my $id ( @order ) {
        $subject_attr_misc_updates{$id}->perform_update;
        push @misc_update_indels, $subject_attr_misc_updates{$id};
    }
    $self->_indels(\@misc_update_indels);

    $self->status_message('Execute INDELS...done');
    return 1;
}

sub _execute_report {
    my $self = shift;

    my $misc_updates = $self->_misc_updates;
    if ( not @$misc_updates ) {
        return print 'RECONCILE MISC UPDATE FOR '.$self->start_from.' TO '.$self->_stop_at."\nNO MISC UPDATES FOUND!\n";
    }

    my (%stats, $msgs);
    my $status = join("\t", (qw/ STATUS SUBJECT_CLASS_NAME SUBJECT_ID SUBJECT_PROPERTY_NAME DESCRIPTION CURRENT_VALUE OLD_VALUE NEW_VALUE /))."\n";
    for my $update ( @{$self->_updates} ) {
        $stats{ATTEMPTED}++;
        $stats{ $update->result }++;
        $status .= $update->status."\n";
        if ( $update->result eq 'FAIL' ) {
            $msgs .= $update->id." '".$update->error_message."'\n";
        }
        elsif ( $update->result eq 'SKIP' ) {
            $msgs .= $update->id." '".$update->status_message."'\n";
        }
    }

    for my $indel ( @{$self->_updates} ) {
        $stats{ATTEMPTED}++;
        $stats{ $indel->result }++;
        $status .= $indel->status."\n";
        if ( $indel->result eq 'FAIL' ) {
            $msgs .= $indel->id." '".$indel->error_message."'\n";
        }
        elsif ( $indel->result eq 'SKIP' ) {
            $msgs .= $indel->id." '".$indel->status_message."'\n";
        }
    }

    return print join(
        "\n\n",
        'RECONCILE MISC UPDATE FOR '.$self->start_from.' TO '.$self->_stop_at,
        "STATS:\n".join("\n", map { sprintf('%-10s => %s', $_, $stats{$_}) } sort keys %stats),
        "STATUS:\n$status",
        "FAILURES and SKIPS [These will remain unreconciled until addressed!]:\n".( $msgs // "NONE :)\n" ),
    );
}

1;

