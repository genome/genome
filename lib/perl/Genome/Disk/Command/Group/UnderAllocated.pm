package Genome::Disk::Command::Group::UnderAllocated;

use strict;
use warnings;

use Genome;

class Genome::Disk::Command::Group::UnderAllocated {
    is => 'Command::V2',
    has_optional => [
        disk_group_names => {
            is => 'Text',
            is_many => 1,
            doc => 'disk groups to be checked',
            default => [ Genome::Config::get('disk_group_alignments'), Genome::Config::get('disk_group_models') ],
        },
        percent_tolerance => {
            is => 'Number',
            default => 1,
            doc => 'How far over the allocation allowed',
        },
    ],
};

sub help_brief {
    return "Finds volumes that are under-allocated (used space exceeds allocated space) and reports them";
}

sub help_synopsis { help_brief() }
sub help_detail { help_brief() }

sub execute {
    my $self = shift;

    my %under_allocated_volumes;
    my %under_allocated_allocations;
    my $tolerance = $self->percent_tolerance;

    # Why yes, I do like my if blocks and for loops nested. Thank you for noticing.
    for my $group ($self->disk_group_names) {
        my @volumes = Genome::Disk::Volume->get(disk_group_names => $group, disk_status => 'active', can_allocate => 1);
        next unless @volumes;

        for my $volume (@volumes) {
            my $allocated = $volume->allocated_kb;
            my $used = $volume->used_kb;
            my $percent_allocated = $volume->percent_allocated;
            my $percent_used = $volume->percent_used;
            my $difference = $percent_used - $percent_allocated;

            if ($used > $allocated and $difference > $tolerance) {
                push @{$under_allocated_volumes{$group}},
                    "Volume " . $volume->mount_path .
                    " using $used kB ($percent_used \%) but only " .
                    "$allocated kB ($percent_allocated \%) allocated";
            }
        }
    }

    my $report = $self->_create_report(\%under_allocated_volumes);
    if ($report) {
        $self->status_message($report);
        return;
    }

    return 1;
}

sub _create_report {
    my ($self, $under_allocated_volumes) = @_;

    my $report;
    for my $group (sort keys %{$under_allocated_volumes}) {
        $report .= "Group $group\n";
        for my $volume (@{$under_allocated_volumes->{$group}}) {
            $report .= "\t$volume\n";
        }
        $report .= "\n";
    }
    return $report;
}

1;

