package Genome::Disk::Command::Allocation::ArchiveExpired;

use strict;
use warnings;
use Genome;
use Carp;

class Genome::Disk::Command::Allocation::ArchiveExpired {
    is => 'Command::V2',
    has => [
    disk_group_name => {
            is => 'Text',
            doc => 'Name for disk group to search against',
            is_optional=>1,
      },
      archive_time => {
        is=>'Time',
        doc=>'Time before which allocations should be archived',
        default_value => Date::Format::time2str(UR::Context->date_template, time()),
      }
    ],
    doc => 'archive allocations where the archive_after_time is prior to now.',
};

sub help_detail {
    return 'archives allocations where the archive_after_time is prior to now.';
}

sub help_brief {
    return 'archives the expired allocations';
}

sub execute {
    my $self = shift;

    $self->status_message("Gathering allocations...");


    my %search_params;
    $search_params{disk_group_name} = $self->disk_group_name if $self->disk_group_name;

    my @volumes = grep {!$_->is_archive} Genome::Disk::Volume->get();

    my @allocations = Genome::Disk::Allocation->get(%search_params,
                                                mount_path=>[map {$_->mount_path} @volumes],
                                                kilobytes_requested=>{operator=>'>=', value=>1024**2},
                                                archive_after_time=>{operator=>'<', value=>$self->archive_time});

    $self->status_message("Found " . scalar @allocations . " allocations to archive");

    return 1 if (@allocations == 0);

    $self->debug_message("Starting archive command on these allocations.");

    my $archive_cmd = Genome::Disk::Command::Allocation::Archive->create(allocations=>\@allocations);

    my $rv = $archive_cmd->execute();
    unless ($rv) {
        Carp::confess "Could not execute archive allocation command.";
    }

    $self->status_message("Done archiving, exiting...");
    return 1;
}

sub get_time {
  my $now = Date::Format::time2str(UR::Context->date_template, time());
}

1;

