use strict;
use warnings;

package Genome::Utility::Test;
use base 'Test::Builder::Module';

use Exporter 'import';
our @EXPORT_OK = qw(diff_ok sub_test);

use Carp qw(croak);
use File::Compare qw(compare);
use Test::More;

my $tb = __PACKAGE__->builder;

sub sub_test($$) {
    my ($desc, $code) = @_;
    $tb->ok($code->(), $desc);
}

sub diff_ok($$;%) {
    my ($file_1, $file_2, %o) = @_;
    use Data::Dumper;

    my $name = delete $o{name};
    my $filter = delete $o{filter};
    my $diag = delete $o{diag} // 1;
    my $test = delete $o{test} // 1;

    my @k = keys %o;
    if (@k) {
        croak 'unexpected options passed to diff_ok: ' . join(', ', @k);
    }


    my @filters;
    if (defined $filter) {
        if (ref($filter) eq 'ARRAY') {
            @filters = @$filter;
        } else {
            @filters = ($filter);
        }
    }

    my @compare_args = ($file_1, $file_2);
    if (@filters) {
        push @compare_args, sub {
            for my $filter (@filters) {
                map { $_ =~ s/$filter// } @_;
            }
            my $c = ($_[0] ne $_[1]);
            if ($c == 1 && $diag && $test) {
                $tb->diag("First diff:\n--- " . $file_1 . "\n+++ " . $file_2 . "\n- " . $_[0] . "+ " . $_[1]);
            }
            return $c;
        };
    } else {
        push @compare_args, sub {
            my $c = ($_[0] ne $_[1]);
            if ($c == 1 && $diag && $test) {
                $tb->diag("First diff:\n--- " . $file_1 . "\n+++ " . $file_2 . "\n- " . $_[0] . "+ " . $_[1]);
            }
            return $c;
        }
    }

    if ($test) {
        return $tb->ok(compare(@compare_args) == 0, $name);
    } else {
        return (compare(@compare_args) == 0 ? 1 : 0);
    }
}

1;

__END__

=pod

=head1 NAME

Genome::Utility::Test

=head1 SYNOPSIS

    use Genome::Utiltiy::Test qw(diff_ok sub_test);

    sub_test('this diffs something' => sub {
        diff_ok($file_1, $file_1);
    });

=head1 METHODS

=item sub_test

Mimics Test::More's subtest since Ubuntu 10.04, which we run, does not have a
Test::More recent enough to have subtest support.

=item diff_ok

diff_ok use File::Compare with a few conveiences.

diff_ok($file_1, $file_2, name => '', diag => 0, test => 0);

=over4

=head2 OPTIONS

=item name

Specify a name for a test but requires 'name' key.

=item diag

Disable diag output when a diff is encountered. Added this in case people want to surpress any output.

=item test

Disable test usage, just return status. Added this so I could test diff_ok.
