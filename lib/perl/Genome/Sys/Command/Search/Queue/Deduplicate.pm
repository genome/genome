package Genome::Sys::Command::Search::Queue::Deduplicate;

class Genome::Sys::Command::Search::Queue::Deduplicate {
    is => 'Command',
    doc => 'deduplicate search queue',
};

sub help_detail {
    return 'Iterate through search queue and deduplicate based subject class and ID.';
}

sub execute {
    return Genome::Search::Queue->dedup();
}

1;
