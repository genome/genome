package Genome::Qc::Tool::Picard::MarkDuplicates;

use strict;
use warnings;
use Genome;

class Genome::Qc::Tool::Picard::MarkDuplicates {
    is => 'Genome::Qc::Tool::Picard',
};

sub supports_streaming {
    return 1;
}

sub output_file_accessor {
    return 'metrics_file';
}

sub metrics {
    return (
        pct_duplicate_reads => {
            picard_metric => 'PERCENT_DUPLICATION',
        },
        number_of_optical_duplicates => {
            picard_metric => 'READ_PAIR_OPTICAL_DUPLICATES',
        },
        estimated_library_size => {
            picard_metric => 'ESTIMATED_LIBRARY_SIZE',
        },
    );
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::MarkDuplicates';
}

1;
