package Genome::Utility::Email;

use strict;
use warnings;

use MIME::Lite;
use Genome::Sys::User;

sub send {
    my %params = @_;

    _fill_in_defaults(\%params);

    _check_required_params(\%params);

    _convert_listrefs_to_strings(\%params);

    my $mime_params = _convert_input_params_to_mime_params(\%params);

    my $msg = MIME::Lite->new( %$mime_params );

    _handle_attachments(\%params, $msg);

    $msg->send();
}

sub _convert_input_params_to_mime_params {
    my $params = shift;

    my %mime_params = map { ucfirst($_) => $params->{$_} } qw(from to subject cc);
    $mime_params{Data} = $params->{body};
    return \%mime_params;
}


sub _check_required_params {
    my $params = shift;

    foreach my $required ( qw( to subject body ) ) {
        die "$required is a required parameter to send()" unless exists $params->{$required};
    }
}

sub _convert_listrefs_to_strings {
    my $params = shift;

    foreach my $is_list ( qw( to from ) ) {
        if (ref($params->{$is_list}) eq 'ARRAY') {
            $params->{$is_list} = join(', ', @{ $params->{$is_list} });
        }
    }
}

sub _fill_in_defaults {
    my $params = shift;

    $params->{from} ||= construct_address();
}

sub _handle_attachments {
    my($params, $msg) = @_;

    return unless exists $params->{attachments};

    my $attachments = $params->{attachments};
    unless (ref($attachments) eq 'ARRAY') {
        $attachments = [ $attachments ];
    }


    foreach my $att ( @$attachments ) {
        _validate_attachment($att);
        $msg->attach(
            Type        => $att->{type},
            Path        => $att->{path},
            Filename    => $att->{filename},
            Disposition => $att->{disposition},
        );
    }
}

sub _validate_attachment {
    my $att = shift;

    foreach my $required ( qw( path ) ) {
        die "$required is a required hash key for an attachment" unless exists $att->{$required};
    }

    $att->{type} ||= 'text/plain';
    $att->{disposition} ||= 'attachment';
    $att->{filename} ||= $att->{path};
}

sub construct_address {
    if (defined($ENV{GENOME_USER_EMAIL})) {
        return $ENV{GENOME_USER_EMAIL};
    }
    my $user = shift;
    $user = defined($user) ? $user : $ENV{USER};
    return join('@', $user, $ENV{GENOME_EMAIL_DOMAIN});
}

1;
