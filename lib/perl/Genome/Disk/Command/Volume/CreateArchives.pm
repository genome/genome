package Genome::Disk::Command::Volume::CreateArchives;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Disk::Command::Volume::CreateArchives {
    is => 'Command::V2',
    doc => 'creates archive volumes for regular volumes, if necessary',
};

sub help_detail { return 'creates archive volumes for regular volumes, if necessary' };
sub help_brief { return help_detail() };
sub help_synopsis { return help_detail() };

sub execute {
    my $self = shift;

    $self->status_message("Creating archive volumes for any volumes that lack one!");

    my %report;
    my $archive_group = Genome::Disk::Group->get(disk_group_name => 'info_archive');
    unless ($archive_group) {
        $archive_group = Genome::Disk::Group->create(
            permissions => '775',
            setgid => '1',
            subdirectory => 'info',
            unix_gid => '10006',
            unix_uid => '10102',
            disk_group_name => 'info_archive',
        );
        unless ($archive_group) {
            die "No archive group found and could not create one!";
        }

        push @{$report{groups}}, $archive_group->disk_group_name;
    }

    my @disk_group_names = ($ENV{GENOME_DISK_GROUP_REFERENCES}, $ENV{GENOME_DISK_GROUP_MODELS}, qw/info_alignments/);
    my @groups = Genome::Disk::Group->get(disk_group_name => \@disk_group_names);
    for my $group (@groups) {
        my @volumes = $group->volumes;
        for my $volume (@volumes) {
            my $archive_mount_path = $volume->archive_mount_path;
            unless ($archive_mount_path) {
                if ($volume->is_archive) {
                    die "Wha? Archive volume " . $volume->mount_path . " is in group " . $group->disk_group_name;
                }
                else {
                    die "Could not figure out archive mount path for volume " . $volume->mount_path;
                }
            }

            my $archive_volume = Genome::Disk::Volume->get(mount_path => $archive_mount_path);
            next if $archive_volume;

            $archive_volume = Genome::Disk::Volume->create(
                total_kb => 1_099_511_627_776,        # one P-P-PETABYTE
                hostname => $volume->hostname,
                can_allocate => 1,
                disk_status => 'active',
                mount_path => $archive_mount_path,
                physical_path => $volume->physical_path,
            );
            unless ($archive_volume) {
                die 'Could not create archive volume';
            }

            my $assignment = Genome::Disk::Assignment->create(
                group => $archive_group,
                volume => $archive_volume
            );
            unless ($assignment) {
                die "Could not add archive volume to archive group!";
            }

            push @{$report{volumes}}, $archive_mount_path;
        }
    }

    if (%report) {
        require Data::Dumper;
        $self->status_message("Create the following items:\n" . Data::Dumper::Dumper(\%report));
    }
    else {
        $self->status_message("No new volumes created");
    }
    
    return 1;
}

