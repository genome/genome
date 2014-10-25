package Genome::Disk::Group;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use Memoize qw(memoize);;
use Module::Find qw(findsubmod usesub);

usesub Genome::Disk::Group::Validate;

class Genome::Disk::Group {
    table_name => 'DISK_GROUP',
    id_by => [
        dg_id => { is => 'Number' },
    ],
    has => [
        disk_group_name => { is => 'Text' },
        permissions => { is => 'Number' },
        setgid => { is => 'Number', is_transient => 1, is_optional => 1 }, #transient during transition from "sticky" column
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
        total_kb => {
            calculate => q(
                my $agg_bx = Genome::Disk::Volume->define_boolexpr(
                    disk_group_names => $self->disk_group_name,
                );
                $self->_uncached_aggregate_sum($agg_bx, 'total_kb');
            ),
        },
        allocated_kb => {
            calculate => q(
                my $agg_bx = Genome::Disk::Allocation->define_boolexpr(
                    disk_group_name => $self->disk_group_name,
                );
                $self->_uncached_aggregate_sum($agg_bx, 'kilobytes_requested');
            ),
        },
        percent_allocated => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ return sprintf("%.2f", ( $allocated_kb / $total_kb ) * 100); },
        },
        unallocated_kb => {
            calculate_from => ['total_kb', 'allocated_kb'],
            calculate => q{ return $total_kb - $allocated_kb },
        },
        used_kb => {
            calculate => q{
                my $sum = 0;
                my $vi = Genome::Disk::Volume->create_iterator(
                    disk_group_names => $self->disk_group_name,
                );
                while (my $v = $vi->next) {
                    $sum += $v->used_kb;
                }
                return $sum;
            },
        },
        percent_used => {
            calculate_from => ['total_kb', 'used_kb'],
            calculate => q{ return sprintf("%.2f", ( $used_kb / $total_kb ) * 100); },
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
    doc => "Represents a disk group (eg, $ENV{GENOME_DISK_GROUP_DEV}), which contains any number of disk volumes",
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

sub validate {
    my $self = shift;
    my @validators = findsubmod Genome::Disk::Group::Validate;
    for (@validators) {
        $_->validate($self);
    }
}

memoize('_resolve_user_name');
sub _resolve_user_name {
    my($self, $uid) = @_;
    return getpwuid($uid);
}

memoize('_resolve_group_name');
sub _resolve_group_name {
    my($self, $gid) = @_;
    return getgrgid($gid);
}

sub _uncached_aggregate_sum {
    my ($self, $agg_bx, $agg_field) = @_;
    my $set = $agg_bx->subject_class_name->define_set($agg_bx);
    $set->__invalidate_cache__("sum($agg_field)");
    return ($set->sum($agg_field) || 0);
}

1;
