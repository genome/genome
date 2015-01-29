package Genome::Notify::Notice;

use strict;
use warnings;

use Genome;

class Genome::Notify::Notice {
    is => ['Genome::Utility::ObjectWithTimestamps','Genome::Utility::ObjectWithCreatedBy'],
    table_name => 'notify.notice',
    has => [
        id => {
            is => 'Text',
        },
        type => {
            is => 'Genome::Notify::NoticeType',
            id_by => 'type_id',
            doc => 'The type of notice',
        },
        header => {
            is => 'Text',
            doc => 'Header (or "subject line") of the notice',
        },
        body => {
            is => 'Text',
            doc => 'Body of the notice',
        },
        subject => {
            is => 'UR::Object',
            id_by => 'subject_id',
            id_class_by => 'subject_class_name',
            doc => 'Entity about which this notice was issued',
        },
        acknowledged => {
            is => 'Boolean',
            doc => 'Whether this notice has been seen or noticed',
        },
        target => {
            is => 'Genome::Sys::User',
            id_by => 'target_id',
            doc => 'The person who is to be notified',
        },
    ],
};

sub issue {
    my $class = shift;
    my @params = @_;

    my $needs_notification;
    my $notice;
    if(my $existing_notice = $class->get(@params)) {
       $notice = $existing_notice;
       $needs_notification = $notice->acknowledged;
       $notice->acknowledged(0);
    } else {
       $needs_notification = 1;
       $notice = $class->create(@params);
    }

    return $notice if $needs_notification;
    return;
}

1;
