use strict;
use warnings;

package Genome::Utility::Test;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(compare_ok sub_test run_ok);

use Carp qw(croak);
use File::Compare qw(compare);
use Test::More;

my $tb = __PACKAGE__->builder;

sub sub_test($$) {
    my ($desc, $code) = @_;
    $tb->ok($code->(), $desc);
}

sub _compare_ok_parse_args {
    # First two args are always the files.
    # Then accept one scalar argument, the name.
    # And one HASH ref, the options.
    # Then validate the option names.
    my $file_1 = shift;
    my $file_2 = shift;

    my $name = (@_ % 2 == 1) ? shift : undef;
    my %o = @_;

    if ($name && $o{name}) {
        die 'duplicate name argument not expected';
    } elsif ($name && !$o{name}) {
        $o{name} = $name;
    }

    my %vo; # validated_o
    $vo{name}    = delete $o{name};
    $vo{filters} = delete $o{filters};
    $vo{diag}    = delete $o{diag} // 1;
    my @k = keys %o;
    if (@k) {
        croak 'unexpected options passed to compare_ok: ' . join(', ', @k);
    }

    my $filters_ref = ref($vo{filters});
    if (defined $vo{filters} && (!$filters_ref || $filters_ref ne 'ARRAY')) {
        $vo{filters} = [$vo{filters}];
    }

    return ($file_1, $file_2, %vo);
}

sub compare_ok {
    my ($file_1, $file_2, %o) = _compare_ok_parse_args(@_);

    my @compare_args = (
        $file_1,
        $file_2,
        sub {
            for my $filter (@{$o{filters}}) {
                map { $_ =~ s/$filter//g } @_;
            }
            my $c = ($_[0] ne $_[1]);
            if ($c == 1 && $o{diag}) {
                $tb->diag("First diff:\n--- " . $file_1 . "\n+++ " . $file_2 . "\n- " . $_[0] . "+ " . $_[1]);
            }
            return $c;
        }
    );

    return $tb->ok(compare(@compare_args) == 0, $o{name});
}

sub run_ok {
    my ($command, $test_name) = @_;

    my $tb = __PACKAGE__->builder;

    my @command = ref $command ? @$command : $command;
    $test_name //= @command > 1 ? join(' ', @command) : $command[0];

    my $exit_zero = (system(@command) == 0);
    $tb->ok($exit_zero, $test_name);

    return $exit_zero;
}

1;

__END__

=pod

=head1 NAME

Genome::Utility::Test

=head1 SYNOPSIS

    use Genome::Utiltiy::Test qw(compare_ok sub_test);

    sub_test('this diffs something' => sub {
        compare_ok($file_1, $file_1);
    });

=head1 METHODS

=item sub_test

Mimics Test::More's subtest since Ubuntu 10.04, which we run, does not have a
Test::More recent enough to have subtest support.

=item compare_ok

compare_ok use File::Compare with a few conveniences.

compare_ok($file_1, $file_2, name => '', diag => 0, test => 0);

=over4

=head2 OPTIONS

=item name

Specify a name for a test but requires 'name' key.

=item diag

Disable diag output when a diff is encountered. Added this in case people want to surpress any output.

=item test

Disable test usage, just return status. Added this so I could test compare_ok.
