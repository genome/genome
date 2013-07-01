use strict;
use warnings;

package Genome::Utility::NamedArgs;

use Exporter 'import';
our @EXPORT_OK = qw(named_args);

use Carp qw(croak);

sub named_args {
    my %argv = @_;

    my $required = delete $argv{required};
    my $optional = delete $argv{optional};
    my $args     = delete $argv{args};

    unless ($required || $optional) {
        croak 'missing argument specification';
    }

    unless ($args) {
        croak 'missing args';
    }

    my @unexpected_argv = keys %argv;
    if (@unexpected_argv) {
        croak 'unexpected named argument(s): ' . join(', ', @unexpected_argv);
    }

    my %o_args = @$args;
    my %o_argv;

    my @missing_required;
    for my $r (@$required) {
        if (exists $o_args{$r}) {
            $o_argv{$r} = delete $o_args{$r};
        } else {
            push @missing_required, $r;
        }
    }
    if (@missing_required) {
        croak 'missing required argument(s): ' . join(', ', @missing_required);
    }

    for my $o (@$optional) {
        $o_argv{$o} = delete $o_args{$o};
    }

    my @unexpected_o_args = keys %o_args;
    if (@unexpected_o_args) {
        croak 'unexpected named argument(s): ' . join(', ', @unexpected_o_args);
    }

    return %o_argv;
}

1;
