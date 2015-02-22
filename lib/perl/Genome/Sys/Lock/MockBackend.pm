package Genome::Sys::Lock::MockBackend;

use strict;
use warnings;

use Test::MockObject qw();

sub new {
    my $self = Test::MockObject->new();
    $self->mock('locks',
        sub {
            my $self = shift;
            $self->{locks} = shift;
            return $self;
        },
    );
    $self->mock('lock',
        sub {
            my $self = shift;
            my %args = @_;
            my $resource_lock = $args{resource_lock};
            if ($self->{locks}{$resource_lock}) {
                return;
            } else {
                $self->{locks}{$resource_lock}++;
                return $resource_lock;
            }
        },
    );
    $self->mock('unlock',
        sub {
            my $self = shift;
            my %args = @_;
            my $resource_lock = $args{resource_lock};
            if ($self->{locks}{$resource_lock}) {
                $self->{locks}{$resource_lock}--;
                return 1;
            } else {
                return;
            }
        },
    );
    $self->mock('has_lock',
        sub {
            my $self = shift;
            my $lock = shift;
            return (defined $self->{locks}{$lock} && $self->{locks}{$lock} > 0);
        },
    );
    $self->mock('release_all',
        sub {
            my $self = shift;
            for my $lock (keys %{$self->{locks}}) {
                $self->unlock(resource_lock => $lock);
            }
        },
    );
    $self->mock('translate_lock_args',   sub { shift; return @_ });
    $self->mock('translate_unlock_args', sub { shift; return @_ });
    $self->set_true('is_mandatory');
    return $self;
}

1;
