package Genome::Site::TGI::Synchronize::ReconcileMiscUpdate;

use strict;
use warnings;

use Genome;

use Date::Format;
use Date::Parse;

class Genome::Site::TGI::Synchronize::ReconcileMiscUpdate {
    is => 'Command::V2',
    has_optional => [
        date => {
            is => 'Text',
            doc => 'Look for updates on this date. Default is to look for updates on the previous day. Format is day-month-year (Date::Format => %d-%b-%y). January 1st 2000 would be "01-JAN-00".',
        },
    ],
    has_optional_transient => [
        _misc_updates => { is => 'Array', default_value => [], },
        stats => { is => 'Hash', default_value => {}, },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return if @errors;

    my $date = $self->date;
    my $time;
    if ( $date ) {
        $time = Date::Parse::str2time($date);
        if ( not $time ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ date /],
                desc => 'Invalid date format => '.$date,
            );
        }
    }
    else {
        $time = time() - 86400;
    }

    $self->date( Date::Format::time2str("%Y-%m-%d", $time) );

    return @errors;
}

sub execute {
    my $self = shift;

    my $execute_updates = $self->_execute_updates;
    return if not $execute_updates;

    my $execute_indels = $self->_execute_indels;
    return if not $execute_indels;

    return $self->_execute_report;
}

sub _execute_updates {
    my $self = shift;

    my @misc_updates = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->get(
        'edit_date like' => $self->date.' %',
        is_reconciled => 0,
        description => 'UPDATE',
        '-order_by' => 'edit_date',
    );
    push @{$self->_misc_updates}, @misc_updates;

    for my $misc_update ( @misc_updates ) {
        $misc_update->perform_update;
    }

    return 1;
}

sub _execute_indels {
    my $self = shift;

    my @misc_updates = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->get(
        'edit_date like' => $self->date.' %',
        is_reconciled => 0,
        'description in'=> [qw/ INSERT DELETE /],
        '-order_by' => 'edit_date',
    );
    push @{$self->_misc_updates}, @misc_updates;

    my @multi_misc_updates;
    foreach my $misc_update ( @misc_updates ) {
        my %multi_misc_update_params = map { $_ => $misc_update->$_ } (qw/ subject_class_name subject_id edit_date description /);
        my $multi_misc_update = Genome::Site::TGI::Synchronize::Classes::MultiMiscUpdate->get(%multi_misc_update_params);
        if ( not $multi_misc_update ) {
            $multi_misc_update = Genome::Site::TGI::Synchronize::Classes::MultiMiscUpdate->create(
                %multi_misc_update_params,
            );
            push @multi_misc_updates, $multi_misc_update;
        }
        $multi_misc_update->add_misc_update($misc_update);
    }

    for my $multi_misc_update ( @multi_misc_updates ) {
        $multi_misc_update->perform_update;
    }

    return 1;
}

sub _execute_report {
    my $self = shift;

    my $misc_updates = $self->_misc_updates;
    if ( not @$misc_updates ) {
        return print 'RECONCILE MISC UPDATE FOR '.$self->date."\nNO MISC UPDATES FOUND!\n";
    }

    my (%stats, $errors);
    my $status = join("\t", (qw/ STATUS SUBJECT_CLASS_NAME SUBJECT_ID SUBJECT_PROPERTY_NAME DESCRIPTION CURRENT_VALUE OLD_VALUE NEW_VALUE /))."\n";
    for my $misc_update ( @$misc_updates ) {
        $stats{ATTEMPTED}++;
        $stats{ $misc_update->result }++;
        $status .= $misc_update->status."\n";
        $errors .= $misc_update->id." '".$misc_update->error_message."'\n" if $misc_update->result eq 'FAILED';
    }

    return print join(
        "\n\n",
        'RECONCILE MISC UPDATE FOR '.$self->date,
        "STATS:\n".join("\n", map { sprintf('%-10s => %s', $_, $stats{$_}) } sort keys %stats),
        "STATUS:\n$status",
        "ERRORS [These will remain unreconciled until addressed!]:\n".( $errors // "NONE :)\n" ),
    );
}

1;

