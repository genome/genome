package Genome::Site::TGI::Synchronize::ReconcileMiscUpdate;

use strict;
use warnings;

use Genome;

use Date::Format;
use Mail::Sender;

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
        _stats => { is => 'Hash', default_value => {}, },
        _status => { is => 'Array', default_value => [], },
        _errors => { is => 'Hash', default_value => {}, },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return if @errors;

    my $date = $self->date;
    if ( $date ) {
        my @tokens = split('-', $date);
        # There is a better way to do this, but it should not be run stand alone very often
        if ( not @tokens or @tokens != 3 or $tokens[0] !~ /^\d{2}$/ or $tokens[1] !~ /^[A-Z]{3}$/ or $tokens[2] !~ /^\d{2}$/ ) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ date /],
                desc => 'Invalid date format => '.$date,
            );
        }
    }
    else {
        my $format = "%d-%b-%y";
        my $now = time();
        my $date = uc(Date::Format::time2str($format, ($now - 86400)));
        $self->date($date);
    }

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Reconcile Misc Update...');

    my $execute_updates = $self->_execute_updates;
    return if not $execute_updates;

    my $execute_indels = $self->_execute_indels;
    return if not $execute_indels;

    $self->_execute_report;
    $self->status_message('Done.');
    return 1;
}

sub _execute_updates {
    my $self = shift;

    $self->status_message('Executing UPDATES...');
    $self->status_message('Date is: '.$self->date);
    my @misc_updates = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->get(
        'edit_date like' => $self->date.' %',
        is_reconciled => 0,
        description => 'UPDATE',
        '-order_by' => 'edit_date',
    );
    push @{$self->_misc_updates}, @misc_updates;
    $self->status_message('Updates found: '.@misc_updates);

    for my $misc_update ( @misc_updates ) {
        $self->status_message( $misc_update->__display_name__ );
        $misc_update->perform_update;
    }

    $self->status_message('UPDATES done');
    return 1;
}

sub _execute_indels {
    my $self = shift;
    $self->status_message('Executing INSERTS/DELETES...');

    $self->status_message('Date is: '.$self->date);
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
    $self->status_message('Inserts/Deletes found: '.@misc_updates);

    for my $multi_misc_update ( @multi_misc_updates ) {
        $self->status_message( $multi_misc_update->__display_name__ );
        $multi_misc_update->perform_update;
    }

    $self->status_message('INSERTS/DELETES done');
    return 1;
}

sub _execute_report {
    my $self = shift;

    my (%stats, @status, %errors);
    for my $misc_update ( @{$self->_misc_updates} ) {
        $stats{attempted}++;
        $stats{ $misc_update->result }++;
        push @status, $misc_update->status;
        $errors{ $misc_update->id } = $misc_update->error_message if $misc_update->result eq 'FAILED';
    }

#    my $sender = Mail::Sender->new({ smtp => 'gscsmtp.wustl.edu', from => 'jlolofie@genome.wustl.edu', });
#    $sender->MailMsg({
#        to => 'apipe-builder@genome.wustl.edu',
#        subject => 'lims -> apipe sync',
#        msg => join("\n",@msg)
#    });

    print Data::Dumper::Dumper( \%stats, \@status, \%errors);

    return 1;
}

1;

