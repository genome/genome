package Genome::Notify::Preference;

use strict;
use warnings;

use Params::Validate qw();

use Genome;

class Genome::Notify::Preference {
    is => ['Genome::Utility::ObjectWithTimestamps'],
    table_name => 'notify.preference',
    has => [
        id => {
            is => 'Text',
        },
        type => {
            is => 'Genome::Notify::NoticeType',
            id_by => 'type_id',
        },
        notifier_class => {
            is_optional => 1,
            is => 'Text',
            doc => 'An override for the default notifier for this type',
        },
        target => {
            is => 'Genome::Sys::User',
            id_by => 'target_id',
            doc => 'The user whose preference this is',
        },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    my $notifier = $self->notifier_class;
    if($notifier and not UNIVERSAL::isa($notifier, 'Genome::Notify::Notifier')) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['notifier_class'],
            desc => 'default_notifier_class must be a Genome::Notify::Notifier',
        );
    }

    return @errors;
}

sub notifier_for_target {
    my $class = shift;
    my %options = Params::Validate::validate(@_, {
        target => { isa => 'Genome::Sys::User' },
        type   => { isa => 'Genome::Notify::NoticeType' },
    });

    my $preference = Genome::Notify::Preference->get(%options);

    return $preference->notifier_class if $preference;
    return $options{type}->default_notifier_class;
}

1;
