use strict;
use warnings;
use above 'Genome';
use Test::More;

no_auto_truncate_message_does_not_truncate_long_text();
auto_truncate_message_truncates_long_text();

done_testing();

sub no_auto_truncate_message_does_not_truncate_long_text {
    my $max_length = Genome::Site::TGI::PSEErrorLogEntry->__meta__->property('message')->data_length;
    my $message = "A"x($max_length + 1);

    my $log = Genome::Site::TGI::PSEErrorLogEntry->create(
        auto_truncate_message => 0,
        message => $message,
    );
    isa_ok($log, 'Genome::Site::TGI::PSEErrorLogEntry', 'log');
    is(length($log->message), ($max_length + 1), "message was not truncated to $max_length");
}
sub auto_truncate_message_truncates_long_text {
    my $max_length = Genome::Site::TGI::PSEErrorLogEntry->__meta__->property('message')->data_length;
    my $message = "A"x($max_length + 1);

    my $log = Genome::Site::TGI::PSEErrorLogEntry->create(
        message => $message,
        auto_truncate_message => 1,
    );
    isa_ok($log, 'Genome::Site::TGI::PSEErrorLogEntry', 'log');
    is(length($log->message), $max_length, "message was truncated to $max_length");
}
