package Genome::Model::Build::Command::Diff;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::Diff {
    is => 'Command::V2',
    has => [
        blessed_build => {
            is => 'Genome::Model::Build',
            doc => "The build which is known to have the correct output.",
        },
        new_build => {
            is => 'Genome::Model::Build',
            doc => "The build that you'd like to know if it has the correct output.",
        },
        _diffs => {
            is => 'HASH',
            doc => "Output from Genome::Model::Build::compare_output function.",
            is_optional => 1,
        },
    ],
};

sub help_brief {
    return 'Determine there are differences in the build directories of two builds.';
}

sub help_detail {
    return <<EOS
    This tool will use Genome::Model::Build::compare_output.
EOS
}


sub execute {
    my $self = shift;

    my $blessed_build = $self->blessed_build;
    my $new_build = $self->new_build;

    $self->status_message(sprintf("Comparing new build %s to blessed build %s...", $new_build->id, $blessed_build->id));
    my %diffs = $blessed_build->compare_output($new_build->id);

    $self->_diffs(\%diffs);

    unless (%diffs) {
        $self->status_message("All files and metrics diffed cleanly!");
    } else {
        $self->status_message($self->diffs_message());
    }

    return 1;
}

sub diffs_message {
    my $self = shift;
    my $diffs = $self->_diffs;
    my $diff_string = "DIFFERENCES FOUND:\n";
    for my $file (sort keys %$diffs) {
        my $reason = $diffs->{$file};
        $diff_string .= "  File: $file\n  Reason: $reason\n";
    }
    return $diff_string;
}
