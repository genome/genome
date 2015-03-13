package Genome::Qc::Tool::BamFile;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::BamFile {
    is => 'Genome::Qc::Tool',
    has => {
        input_file => {
        },
    },
};

sub supports_streaming {
    return 1;
}

sub cmd_line {
    my $self = shift;
    return (qw(samtools view), $self->input_file);
}

1;

