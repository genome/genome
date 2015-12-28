package Genome::Role::Comparable::Command::Diff;

use strict;
use warnings;
use Genome;
use UR::Role;

role Genome::Role::Comparable::Command::Diff {
    has => [
        _diffs => {
            is => 'HASH',
            doc => "Output from Genome::Model::Build::compare_output function.",
            is_optional => 1,
        },
    ],
    requires => ['blessed_object', 'new_object'],
};

sub help_brief {
    return 'Determine if there are differences in the outputs of two comparable objects.';
}

sub help_detail {
    return <<EOS
    This tool will use the compare_output method of the comparable objects.
EOS
}

sub execute {
    my $self = shift;

    my $blessed_object = $self->blessed_object;
    my $new_object = $self->new_object;

    $self->status_message(sprintf("Comparing new object %s to blessed object %s...", $new_object->id, $blessed_object->id));
    my %diffs = $blessed_object->compare_output($new_object->id);

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
    $diff_string .= sprintf(" Comparing new object %s to blessed object %s\n", $self->new_object->id, $self->blessed_object->id);
    for my $file (sort keys %$diffs) {
        my $reason = $diffs->{$file};
        $diff_string .= "  File: $file\n  Reason: $reason\n";
    }
    return $diff_string;
}

sub has_diffs {
    my $self = shift;
    return scalar(keys %{$self->_diffs});
}
1;

