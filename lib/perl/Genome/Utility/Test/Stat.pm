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
    my ($mode, $bit) = (shift, shift);
    my $tb = __PACKAGE__->builder;
    return $tb->ok(_has_bit($mode, $bit), @_);
}

sub hasnt_bit {
    my ($mode, $bit) = (shift, shift);
    my $tb = __PACKAGE__->builder;
    return $tb->ok(!_has_bit($mode, $bit), @_);
}

1;
