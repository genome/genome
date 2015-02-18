package Genome::Notify::NoticeType;

use strict;
use warnings;

use Genome;

class Genome::Notify::NoticeType {
    is => ['UR::Object'],
    table_name => 'notify.notice_type',
    has => [
        id => {
            is => 'Text',
        },
        name => {
            is => 'Text',
            doc => 'The type of notice',
        },
        default_notifier_class => {
            is => 'Text',
            doc => 'The default way to notify a target for notices of this type',
        },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    my $notifier = $self->default_notifier_class;
    if($notifier and not UNIVERSAL::isa($notifier, 'Genome::Notify::Notifier')) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['default_notifier_class'],
            desc => 'default_notifier_class must be a Genome::Notify::Notifier',
        );
    }

    return @errors;
}

1;
