package Genome::Site::TGI::Observers;

use strict;
use warnings;

use Log::Log4perl qw(:easy);

UR::Object::Type->add_observer(
    aspect => 'load',
    callback => sub {
        my $meta = shift;
        my $class_name = $meta->class_name;
        if ($class_name eq 'Genome::ModelGroup') {
            require Genome::Site::TGI::Observers::ModelGroup;
        } elsif ($class_name eq 'Genome::Project') {
            require Genome::Site::TGI::Observers::Project;
        } elsif ($class_name eq 'Command::V1') {
            require Genome::Site::TGI::Observers::Command;
        }
        die $@ if $@;
    },
);

my $log4perl;
UR::Object->add_observer(
    aspect => 'error_message',
    callback => sub {
        my($self, $type, $message) = @_;

        unless ($log4perl) {
            Log::Log4perl->easy_init($ERROR);
            $log4perl = get_logger();
            $log4perl->level($ERROR);

            my $appender = Log::Log4perl::Appender->new(
                "Log::Dispatch::Syslog",
                ident => 'GMS',
                facility => 'syslog',
            );
            $log4perl->add_appender($appender);
        }
        $log4perl->error($message);
    }
);

1;

