package Genome::Disk::Command::Allocation::DoNotArchiveReport;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::DoNotArchiveReport {
    is => 'Command::V2',
};

sub help_detail {
    return 'prints information about allocations that have been ' . 
        'set as do-not-archive grouped by user name';
}
sub help_brief { return 'print information about do-not-archive allocations' };
sub help_synopsis { help_brief() };

sub execute {
    my $self = shift;

    my @allocations = Genome::Disk::Allocation->get(archivable => 0);
    my %report;
    for my $allocation (@allocations) {
        my $user;
        my @notes = $allocation->notes(header_text => 'set to unarchivable');
        if (@notes) {
            $user = $notes[-1]->editor_id;
        }
        else {
            $user = 'unknown';
        }

        $report{$user}{allocations}++;

        my $current_total_kb = $report{$user}{total_kb} || 0;
        $current_total_kb += $allocation->kilobytes_requested;
        $report{$user}{total_kb} = $current_total_kb;
    }

    my @users_by_total_kb = sort { $report{$b}{total_kb} <=> $report{$a}{total_kb} } sort keys %report;
    for my $user (@users_by_total_kb) {
        print join(',', $user, $report{$user}{total_kb}, $report{$user}{allocations}) . "\n";
    }

    return 1;
}

1;

