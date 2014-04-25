package Genome::Site::TGI::Extension::Logger;

use strict;
use warnings;

use Genome::Logger;
use Log::Dispatch::Syslog qw();

my $logger = \&Genome::Logger::logger;
no warnings 'redefine';
*Genome::Logger::logger = sub {
    my $logger = $logger->(@_);
    my $syslog = Log::Dispatch::Syslog->new(
        name => 'syslog',
        min_level => 'debug',
    );
    $logger->add($syslog);
    return $logger;
};

1;
