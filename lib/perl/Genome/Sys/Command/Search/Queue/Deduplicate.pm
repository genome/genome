package Genome::Sys::Command::Search::Queue::Deduplicate;

class Genome::Sys::Command::Search::Queue::Deduplicate {
    is => 'Command',
    doc => 'deduplicate search queue',
};

sub help_detail {
    return 'Iterate through search queue and deduplicate based subject class and ID.';
}

sub execute {
    my $self = shift;

    my %seen;
    my $qi = Genome::Search::Queue->create_iterator();
    while (my $q = $qi->next) {
        my $s_class = $q->subject_class;
        my $s_id = $q->subject_id;
        if ($seen{$s_class}{$s_id}) {
            $q->delete;
        } else {
            $seen{$s_class}{$s_id}++;
        }
    }

    return 1;
}

1;
