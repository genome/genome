package Genome::Report::Command::Email;

use strict;
use warnings;

use Genome;

class Genome::Report::Command::Email {
    is => 'Genome::Report::Command',
    has => [ 
        xsl_files => {
            is => 'Text',
            doc => 'Xslt file(s) to use to transform the report - separate by commas.',
        },
        to => {
            is => 'Text',
            doc => 'Report recipient(s) - separate by commas.',
        },
    ],
    has_optional => [
        from => {
            is => 'Text',
            doc => "Sender of the email.  Defaults to YOU\@$ENV{GENOME_EMAIL_DOMAIN}.",
        },
        replyto => {
            is => 'Text',
            doc => "Reply to for email.  Defaults to YOU\@$ENV{GENOME_EMAIL_DOMAIN}.",
        },
    ],
};

#< Helps >#
sub help_detail {
    return <<EOS;
    Transforms a report with the xsl files, then emails it.
EOS
}

#< Command >#
sub execute {
    my $self = shift;

    # Report
    unless ( $self->report ) {
        $self->error_message("Can't get report.  See above error.");
        $self->delete;
        return;
    }

    my $confirmation = Genome::Report::Email->send_report(
        report => $self->report,
        xsl_files => [ split(',', $self->xsl_files) ],
        to => $self->to,
        from => $self->from,
        replyto => $self->replyto,
    );

    unless ( $confirmation ) {
        $self->error_message("Can't email report.");
        return;
    }

    $self->status_message("Sent report.");

    return 1;
}

1;

#$HeadURL$
#$Id$
