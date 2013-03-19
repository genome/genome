use strict;
use warnings;

package Genome::Utility::Test;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(compare_ok sub_test run_ok capture_ok);

use Carp qw(croak);
use File::Compare qw(compare);
use IPC::System::Simple qw(capture);
use Test::More;
use File::Spec qw();


my %ERRORS = (
    'REPLACE_ARRAY_REF' => q('replace' value should be an ARRAY ref),
);
sub ERRORS {
    my ($class, $key) = @_;
    return $ERRORS{$key};
}

sub sub_test($$) {
    my ($desc, $code) = @_;
    my $tb = __PACKAGE__->builder;
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
    $vo{replace} = delete $o{replace};
    $vo{diag}    = delete $o{diag} // 1;
    my @k = keys %o;
    if (@k) {
        croak 'unexpected options passed to compare_ok: ' . join(', ', @k);
    }

    my $filters_ref = ref($vo{filters});
    if (defined $vo{filters} && (!$filters_ref || $filters_ref ne 'ARRAY')) {
        $vo{filters} = [$vo{filters}];
    }

    my $replace_ref = ref($vo{replace});
    if (defined $vo{replace} && (!$replace_ref || $replace_ref ne 'ARRAY')) {
        die sprintf(q(%s: %s\n), $ERRORS{REPLACE_ARRAY_REF}, $vo{replace});
    }

    return ($file_1, $file_2, %vo);
}

sub compare_ok {
    my ($file_1, $file_2, %o) = _compare_ok_parse_args(@_);

    my $tb = __PACKAGE__->builder;

    my @compare_args = (
        $file_1,
        $file_2,
        sub {
            for my $pr (@{$o{replace}}) {
                my ($pattern, $replacement) = @{$pr};
                map { $_ =~ s/$pattern/$replacement/ } @_;
            }
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

sub capture_ok {
    my ($command, $test_name) = @_;

    my $tb = __PACKAGE__->builder;

    my @command = ref $command ? @$command : $command;
    $test_name //= @command > 1 ? join(' ', @command) : $command[0];

    my @output = eval { capture(@command) };
    my $error = $@;
    my $exit_zero = ($? == 0);
    $tb->ok($exit_zero, $test_name) or diag $error, @output;

    if (wantarray) {
        return ($exit_zero, @output);
    } else {
        return $exit_zero;
    }
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

sub data_dir_ok {
    my ($class, $package) = @_;
    (my $dirname = $package) =~ s/::/-/g;
    my $dirpath = File::Spec->join($ENV{GENOME_TEST_INPUTS}, $dirname);
    my $tb = __PACKAGE__->builder;
    $tb->ok(-d $dirpath, "data_dir exists: $dirpath");
    return $dirpath;
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
