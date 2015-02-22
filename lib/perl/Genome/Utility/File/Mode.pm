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

# as a function it creates an object and as instance method it returns mode
sub mode {
    my $path = shift;

    if (ref $path) {
        return $path->{mode};
    }

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
    return ($mode & $bit) == $bit;
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

sub path { shift->{path} }
sub perm { shift->mode & 07777 }

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

my %bit = (
    setuid           => S_ISUID,
    setgid           => S_ISGID,
    other_rwx        => S_IRWXO,
    other_readable   => S_IROTH,
    other_writable   => S_IWOTH,
    other_executable => S_IXOTH,
    group_rwx        => S_IRWXG,
    group_readable   => S_IRGRP,
    group_writable   => S_IWGRP,
    group_executable => S_IXGRP,
    user_rwx         => S_IRWXU,
    user_readable    => S_IRUSR,
    user_writable    => S_IWUSR,
    user_executable  => S_IXUSR,
    all_rwx          => S_IRWXU | S_IRWXG | S_IRWXO,
    all_readable     => S_IRUSR | S_IRGRP | S_IROTH,
    all_writable     => S_IWUSR | S_IWGRP | S_IWOTH,
    all_executable   => S_IXUSR | S_IXGRP | S_IXOTH,
);
sub bit { %bit }
for my $name (keys %bit) {
    my $bit = $bit{$name};

    my $is = sub { shift->has_bit($bit) };
    install_sub({
        code => $is,
        into => __PACKAGE__,
        as   => 'is_' . $name
    });

    my $add = sub {
        my $self = shift;
        my $mode = $self->{mode} | $bit;
        return $self->set_mode($mode);
    };
    install_sub({
        code => $add,
        into => __PACKAGE__,
        as   => 'add_' . $name
    });

    my $rm = sub {
        my $self = shift;
        my $mode = $self->{mode} & ~$bit;
        return $self->set_mode($mode);
    };
    install_sub({
        code => $rm,
        into => __PACKAGE__,
        as   => 'rm_' . $name
    });
}

sub set_mode {
    my $self = shift;
    my $mode = shift;
    chmod $mode, $self->{path} or croak sprintf('cannot chmod: %s: %s', $self->{path}, $!);
    $self->restat();
    return $self;
}

1;
