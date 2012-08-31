package Genome::Report::Email;
#:adukes xsl and xls are both used here, needs to be fixed, but just a simple typo

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use Mail::Sender;
use Email::Valid;

class Genome::Report::Email {
};

sub send_report {
    my ($class, %params) = @_;

    # XSL Files
    unless ( $params{xsl_files} ) {
        $class->error_message('No XSL files to transform report to email');
        return;
    }
    for my $xsl_file ( @{$params{xsl_files}} ) {
        unless ( Genome::Sys->validate_file_for_reading($xsl_file) ) {
            $class->error_message("Error with xsl file.  See above");
            return;
        }
    }

    # Report
    my $report = delete $params{report};
    unless ( $report ) {
        $class->error_message('Report is required to be able to email it');
        return;
    }
    
    # Email addresses
    my $to = delete $params{to};
    $class->_validate_to_addressees($to)
        or return;
    my $from = ( defined $params{from} ) 
    ? delete $params{from}
    : Genome::Config->user_email;
    $class->_validate_email_address_string('from', $from)
        or return;
    my $reply_to = ( defined $params{replyto} ) 
    ? delete $params{replyto}
    : Genome::Config->user_email;
    $class->_validate_email_address_string('reply to', $reply_to)
        or return;
    
    eval {
        my $sender = Mail::Sender->new({
                smtp => 'gscsmtp.wustl.edu',
                to => $to, 
                from => $from,
                replyto => $reply_to,
                subject => Genome::Utility::Text::capitalize_words( $report->description ),
                multipart => 'related',
                on_errors => 'die',
            }) or die "Can't create Mail::Sender";

        $sender->OpenMultipart;

        for my $xsl_file ( @{$params{xsl_files}} ) {
            # Transform
            my $xslt = Genome::Report::XSLT->transform_report(
                report => $report,
                xslt_file => $xsl_file,
            );
            unless ( $xslt ) {
                $sender->Cancel;
                die $class->error_message("Can't tranform report with xsl file ($xsl_file)");
            }
            # Attach
            $sender->Part({ctype => 'multipart/alternative'});
            $sender->Part({
                    ctype => $xslt->{media_type},
                    disposition => 'NONE',
                    msg => $xslt->{content},
                });
        }

        $sender->EndPart("multipart/alternative");

        if ( $params{attachments} ) {
            for my $attachment ( @{$params{attachments}} ) {
                $sender->Attach($attachment);
            }
        }

        $sender->Close;
    };

    if ( $@ ) {
        $class->error_message("Error configuring email: $@");
        return;
    }

    return 1;
}

#< Email Address >#
sub _validate_email_address_string {
    my ($class, $type, $address_string) = @_;

    unless ( defined $address_string ) {
        $class->error_message("No addresses specified to send report *$type*.");
        return;
    }

    for my $addy ( split(',', $address_string) ) {
        my $valid_addy = eval {
            Email::Valid->address(
                -address => $addy,
                -mxcheck => 1,
            );
        };
        unless ( $valid_addy ) {
            $class->error_message("Error in *$type* email address: $addy");
            return;
        }
    }

    return 1;
}

sub _validate_to_addressees {
    my ($class, $to) = @_;

    unless ( $to ) {
        $class->error_message("To email report, there must be at least one *to* addressees");
        return;
    }

    unless ( ref($to) ) { 
        return $class->_validate_email_address_string('to', $to);
    }

    unless ( @$to ) {
        $class->error_message("To email report, there must be at least one *to* addressees");
        return;
    }
    for my $addressee ( @$to ) {
        $class->_validate_email_address_string('to', $addressee)
            or return;
    }

    return 1;
}

1;

=pod

=head1 Name

Genome::Report::Email

=head1 Synopsis

Email a Genome::Report

=head1 Usage

 use Genome;

 # Get or generate a report...
 my $report = Genome::Report->create_report_from_directory(...);

 # Transform
 my $confirmation = Genome::Report::Email->send_report(
    # Required
    report => $report, # Genome::Report to send
    xslt_files => [...], # transforms report, then adds to email
    to => 'social@genome.wustl.edu', # whom to send the report, separate by commas
    # Optional
    from => 'me@gmail.com', # optional; defaults to current user, separate by commas
    replyto => 'noreply', # optional; defaults to current user, separate by commas
    attachments => [ # images/files to attach
    { # Ex:
        description => 'GC Logo GIF',
        ctype => 'image/jpeg',
        encoding => 'base64',
        disposition => "inline; filename=\"genome_center_logo.gif\";\r\nContent-ID: <footerimg>",
        file => '/gscmnt/839/info/medseq/images/genome_center_logo.gif'
    },
    # ...
    ],
 );

 unless ( $confirmation ) {
    die "Can't send report\n";
 }
 
 ...
 
=head1 Public Methods

=head2 transform_report

 my $string = Genome::Report::XSLT->transform_report(report => $report, xslt_file => $xslt_file);

=over

=item I<Synopsis>   Takes a report and email xslt file (as a hash), and returns the transformed report as a string

=item I<Arguments>  report (Genome::Report), xslt_files (array of xslt files), addressees (to, from, replyto), image files (array of hashes)

=item I<Returns>    boolean

=back

=head1 Disclaimer

Copyright (C) 2009 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@genome.wustl.edu>

=cut


#$HeadURL$
#$Id$
