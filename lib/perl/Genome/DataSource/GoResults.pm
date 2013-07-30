package Genome::DataSource::GoResults;

use strict;
use warnings;
use Genome;

class Genome::DataSource::GoResults {
    is => [ 'UR::DataSource::FileMux', 'Genome::DataSource::FileMuxDirMustExist', 'UR::Singleton' ],
};

sub delimiter { "," }

sub column_order {
    [ 'go_id', 'interpro_result_id', 'start', 'end', 'transcript_name', 'name', 'term_type' ];
}

sub required_for_get {
    [ 'data_directory', 'chrom_name' ];
}

sub file_resolver {
    my ($data_directory, $chrom_name) = @_;
    my $path = $data_directory . "/go_results/chromosome_" . $chrom_name . ".csv";
    return __PACKAGE__->directory_must_exist($path);
}

1;

