package Genome::Command::Search;

use strict;
use warnings;


class Genome::Command::Search {
    is => 'Command::V2',

    has => [
        target => {
            is => 'Text',
            shell_args_position => 1,
            is_many => 1,
            doc => 'The string, id, or partial id to search for.',
        },

        rows => {
            is => 'Integer',
            is_optional => 1,
            default_value => 15,
            doc => 'Maximum number of rows to return',
        },
    ],
};


sub execute {
    my $self = shift;

    for my $doc ($self->get_docs) {
        printf("%s: %s\n", $doc->{display_type}, $doc->{display_title});
    }

    return 1;
}
sub get_docs {
    my $self = shift;

    my $content = $self->get_content;
    return @{$content->{response}{docs}};
}

sub get_content {
    my $self = shift;

    my $response = Genome::Search->search(join(' ', $self->target),
        {rows => $self->rows});
    unless (defined($response)) {
        die "Invalid response from Search";
    }

    $self->validate_response_content($response->content);

    return $response->content;
}


my $SCORE_WARNING_THRESHOLD = 5;
sub validate_response_content {
    my ($self, $content) = @_;

    if ($content->{response}{numFound} <= 0) {
        die "No results found";
    }

    if ($content->{response}{maxScore} < $SCORE_WARNING_THRESHOLD) {
        warn sprintf("Warning: Low search score (%s < %s).",
            $content->{response}{maxScore}, $SCORE_WARNING_THRESHOLD);
    }
}


1;
