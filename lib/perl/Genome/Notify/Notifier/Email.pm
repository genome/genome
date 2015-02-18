package Genome::Notify::Notifier::Email;

use strict;
use warnings;

use Genome;
use Genome::Utility::Email;

class Genome::Notify::Notifier::Email {
    is => 'Genome::Notify::Notifier',
    has_classwide_constant => [
        name => {
            value => 'e-mail',
        },
    ],
    doc => 'The notice will be e-mailed. The e-mail will not be repeated until the notice is acknowledged.',
};

sub notify {
    my $class = shift;
    my $notice = shift;

    Genome::Utility::Email::send(
        from    => $ENV{GENOME_EMAIL_PIPELINE},
        replyto => $ENV{GENOME_EMAIL_NOREPLY},
        to      => $notice->target->email,
        subject => $notice->header,
        body    => $notice->body,
    );
}

1;
