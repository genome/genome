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

sub qc_metrics_file_accessor {
    return 'metrics_file';
}

sub gmt_class {
    return 'Genome::Model::Tools::Picard::MarkDuplicates';
}

sub output_file {
    return Genome::Sys->create_temp_file_path();
}

1;
