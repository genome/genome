package Genome::DataSource::InterproResults;

use strict;
use warnings;
use Genome;

class Genome::DataSource::InterproResults {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton' ],
};

sub delimiter { 
    return "," 
}

sub column_order {
    return [qw(
        interpro_id
        start
        stop
        transcript_name
        rid
        setid
        parent_id
        name 
        interpro_note
    )]
}

sub constant_values { 
    [ 'data_directory', 'chrom_name' ]; 
};

sub sort_order {
    [qw/ start stop interpro_id /];
}

sub required_for_get { 
    [ 'data_directory', 'chrom_name' ]; 
}

sub file_resolver {
    my ($data_directory, $chrom_name) = @_;
    my $path = $data_directory . "/interpro_results/chromosome_" . $chrom_name . ".csv";
    return $path;
}

1;

