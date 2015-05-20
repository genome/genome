package Genome::Disk::Command::Group::AvailableSpace;

use strict;
use warnings;
use Genome;
use Genome::Utility::Email;

class Genome::Disk::Command::Group::AvailableSpace {
    is => 'Command::V2',
    has_optional => [
        disk_group_names => {
            is => 'Text',
            doc => 'comma delimited list of disk groups to be checked',
        },
    ],
};

my %minimum_space_for_group = (
    Genome::Config::get('disk_group_dev') => 512_000,                # 500MB
    Genome::Config::get('disk_group_references') => 1_073_741_824,      # 1TB
    Genome::Config::get('disk_group_alignments') => 12_884_901_888,    # 12TB
    Genome::Config::get('disk_group_models') => 25_769_803_776  # 24TB
);

sub help_brief {
    return "Sums up the unallocated space for every volume of a group";
}

sub help_synopsis {
    help_brief();
}

sub help_detail {
    help_brief();
}

sub execute {
    my $self = shift;

    my @disk_groups;
    if ($self->disk_group_names) {
        @disk_groups = split(',', $self->disk_group_names);
    }
    else {
        @disk_groups = keys %minimum_space_for_group;
    }

    my $group_is_low = 0;
    my @reports;
    for my $group_name (@disk_groups) {
        my $group = Genome::Disk::Group->get(disk_group_name => $group_name);
        next unless $group;

        my @volumes = grep { $_->can_allocate == 1 and $_->disk_status eq 'active' } $group->volumes;
        next unless @volumes;

        my $sum;
        for my $volume (@volumes) {
            my $space = $volume->soft_limit_kb - $volume->allocated_kb;
            $sum += $space unless $space < 0;  # I've learned not to trust the system to be consistent
        }

        my $sum_gb = $self->kb_to_gb($sum);
        my $sum_tb = $self->kb_to_tb($sum);

        if (exists $minimum_space_for_group{$group_name} and $sum < $minimum_space_for_group{$group_name}) {
            my $min = $minimum_space_for_group{$group_name};
            $group_is_low = 1;
            my $min_gb = $self->kb_to_gb($min);
            my $min_tb = $self->kb_to_tb($min);
            push @reports, "Disk group $group_name has $sum KB ($sum_gb GB, $sum_tb TB) of free space, " .
                "which is below the minimum of $min KB ($min_gb GB, $min_tb TB). Either free some disk or request more!";
        }
        else {
            push @reports, "Disk group $group_name has $sum KB ($sum_gb GB, $sum_tb TB) of free space";
        }
    }

    $self->status_message(join("\n", @reports));

    return !$group_is_low;
}

sub kb_to_gb {
    my ($class, $kb) = @_;
    return int($kb / (2**20));
}

sub kb_to_tb {
    my ($class, $kb) = @_;
    return int($kb / (2**30));
}

1;

