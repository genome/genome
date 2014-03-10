package Genome::Site::TGI::Extension::Logger;

use Genome::Logger;
use Log::Dispatch::Syslog qw();

sub import {
    my $logger = Genome::Logger->logger();
    my $syslog = Log::Dispatch::Syslog->new(
        name => 'syslog',
        min_level => 'debug',
    );
    $logger->add($syslog);
}

1;
