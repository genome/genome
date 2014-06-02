use strict;
use warnings;

package Genome::Utility::Test::Stat;
use parent 'Test::Builder::Module';

use Exporter 'import';
use Fcntl ':mode';
use List::MoreUtils qw(any);

our @EXPORT_OK = qw(has_bit hasnt_bit);

my @DATA = (
    [ 0, [S_IRWXO, S_IROTH, S_IWOTH, S_IXOTH, S_ISUID, S_ISGID]], # i = 0, i.e. first bit
    [ 3, [S_IRWXG, S_IRGRP, S_IWGRP, S_IXGRP]],                   # i = 1, i.e. second bit
    [ 6, [S_IRWXU, S_IRUSR, S_IWUSR, S_IXUSR]],                   # i = 2, i.e. third bit
);

my %names = (
    S_IRWXO, 'other %s read, write, execute',
    S_IROTH, 'other %s read',
    S_IWOTH, 'other %s write',
    S_IXOTH, 'other %s execute',
    S_ISUID, '%s setuid',
    S_ISGID, '%s setgid',
    S_IRWXG, 'group %s read, write, execute',
    S_IRGRP, 'group %s read',
    S_IWGRP, 'group %s write',
    S_IXGRP, 'group %s execute',
    S_IRWXU, 'user %s read, write, execute',
    S_IRUSR, 'user %s read',
    S_IWUSR, 'user %s write',
    S_IXUSR, 'user %s execute',
);

sub _find {
    my $bit = shift;
    for (@DATA) {
        my $shift = $_->[0];
        my @bits  = @{$_->[1]};
        if (any { $bit == $_ } @bits) {
            return $shift;
        }
    }
}

sub _has_bit {
    my ($mode, $bit) = (shift, shift);
    return ($mode & $bit) >> _find($bit);
}

sub has_bit {
    my ($mode, $bit, $name) = @_;
    unless ($name) {
        $name = sprintf($names{$bit}, 'has');
    }
    my $tb = __PACKAGE__->builder;
    return $tb->ok(_has_bit($mode, $bit), $name);
}

sub hasnt_bit {
    my ($mode, $bit, $name) = @_;
    unless ($name) {
        $name = sprintf($names{$bit}, 'does not have');
    }
    my $tb = __PACKAGE__->builder;
    return $tb->ok(!_has_bit($mode, $bit), $name);
}

1;
