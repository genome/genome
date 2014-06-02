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
    S_ISUID, 'setuid',
    S_ISGID, 'setgid',
    S_IRWXO, 'other_rwx',
    S_IROTH, 'other_readable',
    S_IWOTH, 'other_writable',
    S_IXOTH, 'other_executable',
    S_IRWXG, 'group_rwx',
    S_IRGRP, 'group_readable',
    S_IWGRP, 'group_writable',
    S_IXGRP, 'group_executable',
    S_IRWXU, 'user_rwx',
    S_IRUSR, 'user_readable',
    S_IWUSR, 'user_writable',
    S_IXUSR, 'user_executable',
);
for my $bit (keys %names) {
    my $name = $names{$bit};

    my $is = sub { shift->has_bit($bit) };
    install_sub({
        code => $is,
        into => __PACKAGE__,
        as   => 'is_' . $name
    });
}

1;
