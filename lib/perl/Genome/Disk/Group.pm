package Genome::Disk::Group;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Disk::Group {
    table_name => 'DISK_GROUP',
    id_by => [
        dg_id => { is => 'Number' },
    ],
    has => [
        disk_group_name => { is => 'Text' },
        permissions => { is => 'Number' },
        setgid => { is => 'Number', is_transient => 1 }, #transient during transition from "sticky" column
        subdirectory => { is => 'Text' },
        unix_uid => { is => 'Number' },
        unix_gid => { is => 'Number' },
        user_name => {
            calculate_from => 'unix_uid',
            calculate => q|
                my ($user_name) = $self->_resolve_user_name($unix_uid);
                return $user_name;
            |,
        },
        group_name => {
            calculate_from => 'unix_gid',
            calculate => q| 
                my ($group_name) = $self->_resolve_group_name($unix_gid);
                return $group_name;
            |,
        },
    ],
    has_many_optional => [
        mount_paths => {
            via => 'volumes',
            to => 'mount_path',
        },
        volumes => {
            is => 'Genome::Disk::Volume',
            via => 'assignments',
            to =>  'volume',
        },
        assignments => {
            is => 'Genome::Disk::Assignment',
            reverse_id_by => 'group',
        },
    ],
    data_source => 'Genome::DataSource::Oltp',
    doc => 'Represents a disk group (eg, info_apipe), which contains any number of disk volumes',
};


sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    if (defined($self)) {
        $self->validate;
    }

    return $self;
}

sub get {
    my $class = shift;
    if (wantarray) {
        return $class->_get_many(@_);
    } else {
        return $class->_get_single(@_);
    }
}

sub _get_many {
    my $class = shift;

    my @objs = $class->SUPER::get(@_);
    for my $obj (@objs) {
        $obj->validate;
    }

    return @objs;
}

sub _get_single {
    my $class = shift;

    my $self = $class->SUPER::get(@_);
    if (defined($self)) {
        $self->validate;
    }

    return $self;
}

# TODO This needs to be removed, site-specific
my %VALID_NAMES = (
    info_apipe => 1,
    info_apipe_ref => 1,
    info_alignments => 1,
    info_genome_models => 1,
    apipe_trash => 1,
    research => 1,
    systems_benchmarking => 1,
);
sub validate {
    my $self = shift;

    unless ($ENV{UR_DBI_NO_COMMIT}) {
        unless ($VALID_NAMES{$self->disk_group_name}) {
            confess sprintf("Disk group name (%s) not allowed.",
                $self->disk_group_name);
        }
    }
}

my %user_name_cache;
my %group_name_cache;

# memoizing frontends for user_name and group_name
sub _resolve_user_name {
    my($self, $uid) = @_;

    unless (exists $user_name_cache{$uid}) {
        ($user_name_cache{$uid}) = getpwuid($uid);
    }
    return $user_name_cache{$uid};
}

sub _resolve_group_name {
    my($self,$gid) = @_;

    unless (exists $group_name_cache{$gid}) {
        ($group_name_cache{$gid}) = getgrgid($gid);
    }
    return $group_name_cache{$gid};
}


1;
