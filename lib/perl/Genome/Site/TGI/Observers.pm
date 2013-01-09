package Genome::Site::TGI::Observers;

use strict;
use warnings;

use Log::Log4perl qw(get_logger :levels);

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

        # this should never happen given recent UR updates
        if (not defined $self) {
            no warnings;
            Carp::confess("self is undef, are you using the latest UR?: @_");
        }
        if (not defined $message) {
            no warnings;
            Carp::confess("message is undef, are you using the latest UR?: @_");
        }

        unless ($log4perl) {
            $log4perl = get_logger();
            $log4perl->level($ERROR);

            my $appender = Log::Log4perl::Appender->new(
                "Log::Dispatch::Syslog",
                ident => "GMS $0",
                facility => 'syslog',
            );
            $log4perl->add_appender($appender);
        }
        my $a = ref($self) ? $self->class . ' id('. $self->__display_name__.')' : $self;
        $log4perl->error($a . ': ' . $message);
    }
);

1;

