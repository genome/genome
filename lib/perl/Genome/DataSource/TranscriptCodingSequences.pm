package Genome::DataSource::TranscriptCodingSequences;

use strict;
use warnings;
use Genome;

class Genome::DataSource::TranscriptCodingSequences {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub delimiter {
    return ",";
}

sub column_order {
    [ 'transcript_id', 'sequence' ];
}

sub constant_values {
    [ 'data_directory' ];
};

sub required_for_get {
    [ 'data_directory' ];
}

sub file_resolver {
    my $data_directory = shift;
    my $path = $data_directory . "/transcript_coding_sequences.csv";
    return $path;
}

1;

