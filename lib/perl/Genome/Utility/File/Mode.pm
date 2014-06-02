package Genome::Utility::File::Mode;

use strict;
use warnings;

use Carp qw(croak);
use Exporter qw(import);
use Fcntl qw(:mode);
use File::stat qw(stat);
use List::MoreUtils qw(any);
use Sub::Install qw(install_sub);

our @EXPORT_OK = qw(mode);

sub mode {
    my $path = shift;
    __PACKAGE__->new($path);
}

my @DATA = (
    [ 0, [S_IRWXO, S_IROTH, S_IWOTH, S_IXOTH, S_ISUID, S_ISGID]], # i = 0, i.e. first bit
    [ 3, [S_IRWXG, S_IRGRP, S_IWGRP, S_IXGRP]],                   # i = 1, i.e. second bit
    [ 6, [S_IRWXU, S_IRUSR, S_IWUSR, S_IXUSR]],                   # i = 2, i.e. third bit
);
sub _bit_shift {
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
    return ($mode & $bit) >> _bit_shift($bit);
}

sub new {
    my $class = shift;
    my $path  = shift;
    my $self = bless {
        path => $path,
    }, $class;
    $self->restat();
    return $self;
}

sub restat {
    my $self = shift;
    my $stat = stat($self->{path}) or croak "cannot stat: $!";
    $self->{mode} = $stat->mode;
}

sub has_bit {
    my $self = shift;
    my $bit = shift;
    return _has_bit($self->{mode}, $bit);
}

my %names = (
    S_ISUID, 'is_setuid',
    S_ISGID, 'is_setgid',
    S_IRWXO, 'is_other_rwx',
    S_IROTH, 'is_other_readable',
    S_IWOTH, 'is_other_writable',
    S_IXOTH, 'is_other_executable',
    S_IRWXG, 'is_group_rwx',
    S_IRGRP, 'is_group_readable',
    S_IWGRP, 'is_group_writable',
    S_IXGRP, 'is_group_executable',
    S_IRWXU, 'is_user_rwx',
    S_IRUSR, 'is_user_readable',
    S_IWUSR, 'is_user_writable',
    S_IXUSR, 'is_user_executable',
);
for my $bit (keys %names) {
    my $name = $names{$bit};
    my $sub = sub {
        return shift->has_bit($bit);
    };
    install_sub({
        code => $sub,
        into => __PACKAGE__,
        as   => $name
    });
}

1;
