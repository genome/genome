package Genome::Command::Copy;

use strict;
use warnings FATAL => 'all';

use Genome;

use Carp qw(croak);
use Try::Tiny qw(try catch);

class Genome::Command::Copy {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        source => {
            is => 'UR::Object',
        },
        changes => {
            is => 'Text',
            is_many => 1,
            doc => 'comma-separated list of changes',
        },
    ],
};

sub help_detail {
    return join('',
        q(Any non-delegated, non-ID properties may be specified with an operator and a value.),
        q( Valid operators are '=', '+=', '-=', and '.='; function is same as in Perl.),
        qq(\n\nFor example:\n\n),
        qq(   --changes "name.=-RT101912,foo=bar"\n\n),
        q( A value of 'undef' may be used to pass a Perl undef as the value.  Either `foo=` or `foo=''` can be used to set the value to an empty string.),
    );
}

sub execute {
    my $self = shift;

    my $tx = UR::Context::Transaction->begin();

    my $error;
    my $copy = try {
        my $copy = $self->source->copy();

        for my $change ($self->changes) {
            my ($key, $op, $value) = $change =~ /^(.+?)(=|\+=|\-=|\.=)(.*)$/;
            unless ($key && $op) {
                die $self->error_message("invalid change: $change");
            }

            if ($value eq 'undef') {
                $value = undef;
            }

            if ($op eq '=') {
                $copy->$key($value);
            }
            elsif ($op eq '+=') {
                $copy->$key($copy->$key + $value);
            }
            elsif ($op eq '-=') {
                $copy->$key($copy->$key - $value);
            }
            elsif ($op eq '.=') {
                $copy->$key($copy->$key . $value);
            }
        }

        unless ($tx->commit) {
            die 'commit failed';
        }

        return $copy;
    }
    catch {
        $tx->rollback();
        $error = $_ || 'unknown error';
        return;
    };

    if ($copy) {
        $self->status_message('Created new %s with ID %s', $copy->class, $copy->id);
        return 1;
    }
    else {
        $self->error_message('Failed to create new %s: %s', $self->source->class, $error);
        return;
    }
}
