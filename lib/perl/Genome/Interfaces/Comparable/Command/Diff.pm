package Genome::Interfaces::Comparable::Command::Diff;

use strict;
use warnings;
use Genome;

class Genome::Interfaces::Comparable::Command::Diff {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        _diffs => {
            is => 'HASH',
            doc => "Output from Genome::Model::Build::compare_output function.",
            is_optional => 1,
        },
    ],
};

sub help_brief {
    'compare two objects'
}

sub help_detail {
    "This command will compare the outputs of 2 objects using the objects' compare_output method. Any differences found will be displayed."
}

sub blessed_object {
    die "Must implement blessed_object";
}

sub new_object {
    die "Must implement new_object";
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

