package Genome::Disk::Detail::Allocation::PermissionScheme;

use strict;
use warnings;

use Genome;

use Carp qw(croak);
use File::stat qw(stat);
use Path::Class::Dir qw();


class Genome::Disk::Detail::Allocation::PermissionScheme {
    table_name => 'disk.allocation_permission_scheme',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
    id_by => [
        id => {
            is => 'Text',
            len => 64,
        },
    ],
    has => [
        uid => {
            is_optional => 1,
            is => 'Number',
            default => -1,
            doc => q(Not yet supported.  Numeric user ID.  Default value, -1, means leave user ID unchanged.),
        },
        gid => {
            is_optional => 1,
            is => 'Number',
            default => -1,
            doc => q(Numeric group ID.  Default value, -1, means leave group ID unchanged.),
        },
        min_mode => {
            is_optional => 1,
            is => 'Number',
            default => 0,
            doc => q(Minimum mode allowed.  Default value, 0, means no minimum.),
        },
        max_mode => {
            is_optional => 1,
            is => 'Number',
            default => oct(7777),
            doc => q(Maximum mode allowed.  Default value, oct(7777), means no maximum.),
        },
        file_min_mask => {
            is_optional => 1,
            is => 'Number',
            default => 0,
            doc => q(Bit mask to use for files' minimum mode (as opposed to directories).  Default value, 0, means directories and files would have the same minimum mode.),
        },
    ],
};


sub create {
    my $classname = shift;
    my $tx = UR::Context::Transaction->begin();
    my $self = $classname->SUPER::create(@_);
    if ($self->uid != -1) {
        $tx->rollback();
        croak 'setting uid is not yet supported';
    } else {
        $tx->commit();
    }
    return $self;
}


sub apply {
    my $self = shift;
    my $path = shift;

    my $rv = 1;
    my $dir = Path::Class::Dir->new($path);
    $dir->recurse(
        callback => sub {
            my $obj = shift;
            unless (chown $self->uid, $self->gid, $obj->stringify) {
                warn qq(chown failed: $!);
                $rv = undef;
            }
            my $mode = $self->mode_for($obj->stringify);
            unless (chmod $mode, $obj->stringify) {
                warn qq(chmod failed: $!);
                $rv = undef;
            }
        },
    );

    return $rv;
}


sub mode_for {
    my $self = shift;
    my $path = shift;

    if (-l $path) {
        return oct(777);
    } elsif (! -l $path && -d $path) {
        return $self->dir_mode($path);
    } elsif (! -l $path && -f $path) {
        return $self->file_mode($path);
    } else {
        die qq(non-existant path? $path);
    }
}


sub umask {
    my $self = shift;
    return ~int($self->max_mode) & oct(7777);
}


sub dir_mode {
    my $self = shift;
    my $path = shift;
    my $mode = stat($path)->mode;
    return ($mode | int($self->min_mode)) & ~int($self->umask);
}


sub file_min_mode {
    my $self = shift;
    return int($self->min_mode) & ~int($self->file_min_mask);
}


sub file_mode {
    my $self = shift;
    my $path = shift;
    my $mode = stat($path)->mode;
    return ($mode | int($self->file_min_mode)) & ~int($self->umask);
}

1;
