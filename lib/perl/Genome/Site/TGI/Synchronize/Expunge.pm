package Genome::Site::TGI::Synchronize::Expunge;

use strict;
use warnings;
use Genome;
use Mail::Sender;
use Carp 'confess';

class Genome::Site::TGI::Synchronize::Expunge {
    is => 'Genome::Command::Base',
    has => [
        report => {
            is => 'Hashref',
            is_input => 1,
            is_optional => 0,
            doc => 'Hashref containing objects of interest from Genome::Site::TGI::Synchronze::UpdateApipeClasses',
        },
    ],
};

sub execute {
    my $self = shift;
    my %report = %{$self->report};
    my %expunge_notifications;

    for my $class (keys %report){
        next unless $class =~ m/Genome::InstrumentData/; #only remove instrument data for now
        next if $class eq 'Genome::InstrumentData::Imported'; #imported instrument data doesn't come from LIMS, so skip it
        print "WOULD DELETE $class FOR IDS:\n";
        print Data::Dumper::Dumper($report{$class}->{missing});
        next;
        my @ids = @{$report{$class}->{missing}} if $report{$class}->{missing};
        my @deleted;
        for my $id (@ids){
            my ($successfully_deleted, %expunge_info) = $self->_remove_expunged_object($class, $id);
            push @deleted, $successfully_deleted;
            $self->_merge_expunge_notifications(\%expunge_notifications, \%expunge_info);
        }
        $report{$class}->{deleted} = \@deleted;
    }

    $self->_notify_expunged_objects_owners(%expunge_notifications);
    
    return 1;
}

sub _remove_expunged_object {
    my $self = shift;
    my $class = shift;
    my $id = shift;
    my $expunge_success;
    my %affected_users;

    my $object = $class->get($id);

    $object->delete;

    return ($id, %affected_users);
}

sub _merge_expunge_notifications {
    my $self = shift;
    my $master_notifications = shift;
    my $new_notifications = shift;

    for my $user_name (keys %$new_notifications){
        my %tmp;
        if($master_notifications->{$user_name}){
            %tmp =  (%{$new_notifications->{$user_name}}, %{$master_notifications->{$user_name}});
        }else{
            %tmp = (%{$new_notifications->{$user_name}});
        }
        $master_notifications->{$user_name} = \%tmp;
    }

    return %$master_notifications;
}

sub _notify_expunged_objects_owners{
    my $self = shift;
    my %expunge_notifications = @_;
    
    for my $user_name (keys %expunge_notifications){
        my ($subject, $msg) = $self->_generate_expunged_objects_email_text(%{$expunge_notifications{$user_name}});
        my $sender = Mail::Sender->new({
                smtp    => 'gscsmtp.wustl.edu',
                from    => 'Apipe <apipe-builder@genome.wustl.edu>'
                });
        $sender->MailMsg( { 
                to      => 'Analysis Pipeline <apipebulk@genome.wustl.edu>, ' . "$user_name".'@genome.wustl.edu', 
                cc      => 'Jim Weible <jweible@genome.wustl.edu>, Thomas Mooney <tmooney@genome.wustl.edu>, Scott Smith <ssmith@genome.wustl.edu>',
                subject => "Expunged Instrument Data: $subject", 
                msg     => "LIMS has expunged instrument data used in some of your models.  Existing builds using this data will be abandoned and the model will be rebuilt.  Please contact APipe if you have any questions regarding this process.\n\n$msg", 
                });
    }
}

sub _generate_expunged_objects_email_text {
    my $self = shift;
    my %expunge_notifications_for_user = @_;
    my @display_names;
    
    my $message_text = "";
    for my $instrument_data_info (keys %expunge_notifications_for_user){
        my ($display_name, $instrument_data_id) = split(" ", $instrument_data_info);
        push @display_names, $display_name;

        $message_text .= "Instrument Data: $display_name ($instrument_data_id)\n";
        my @model_ids = map($_, @{$expunge_notifications_for_user{$instrument_data_info}});
        for my $model_id (@model_ids){
            my $model = Genome::Model->get($model_id);
            $message_text .= join(" ", "\tModel:", $model->name, "($model_id)\n");
        }
        $message_text .= "\n";
    }
    my $subject_text = join(", ", @display_names);
    return ($subject_text, $message_text);
}
1;
