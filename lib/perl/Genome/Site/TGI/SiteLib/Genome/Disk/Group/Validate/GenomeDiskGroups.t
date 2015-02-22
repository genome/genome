#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 3;

use Genome qw();
use Genome::Disk::Group qw();
use List::MoreUtils qw(none);
use List::Util qw(shuffle);
use Test::Fatal qw(exception);

do {
    my $disk_group_name = random_string();
    local $ENV{GENOME_DISK_GROUP_MODELS} = $disk_group_name;
    local $ENV{UR_DBI_NO_COMMIT} = 0;
    my $group = Genome::Disk::Group->__define__(
        disk_group_name => $disk_group_name,
        permissions => 0,
        subdirectory => 'info',
        unix_uid => 0,
        unix_gid => 0,
    );
    ok( ! exception { $DB::single = 1; $group->validate },
        'no exception thrown when disk_group_name is GENOME_DISK_GROUP_MODELS'
    );
};

do {
    my $disk_group_name = random_string();
    local $ENV{GENOME_DISK_GROUP_MODELS} = random_string();
    local $ENV{UR_DBI_NO_COMMIT} = 0;

    ok((none { $_ eq $disk_group_name } Genome::Disk::Group::Validate::GenomeDiskGroups::genome_disk_group_names()),
        'disk_group_name not in genome_disk_group_names');

    my $group = Genome::Disk::Group->__define__(
        disk_group_name => $disk_group_name,
        permissions => 0,
        subdirectory => 'info',
        unix_uid => 0,
        unix_gid => 0,
    );
    like( exception { $DB::single = 1; $group->validate },
        qr/not allowed/,
        'exception thrown when disk_group_name is not GENOME_DISK_GROUP_MODELS'
    );
};

sub random_string {
    my @chars = map { (shuffle 'a'..'z', 'A'..'Z', 0..9)[0] } 1..10;
    return join('', @chars);
}
