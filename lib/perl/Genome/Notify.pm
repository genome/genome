package Genome::Notify;

use strict;
use warnings;

use Params::Validate qw();

use Genome;

class Genome::Notify {
    is => 'UR::Singleton',
};

sub notify {
    my $class = shift;
    my %options = Params::Validate::validate(@_, {
        subject => { isa => 'UR::Object' },
        type    => { isa => 'Genome::Notify::NoticeType' },
        header  => { type => Params::Validate::SCALAR },
        body    => { type => Params::Validate::SCALAR },
        target  => { isa => 'Genome::Sys::User' },
    });

    my $notice = Genome::Notify::Notice->issue(%options);
    return unless $notice;

    my $notifier = Genome::Notify::Preference->notifier_for_target(
        type => $options{type},
        target => $options{target},
    );

    $notifier->notify($notice);
}

1;
