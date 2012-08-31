package Genome::DataSource::Transcripts;

use Genome;

class Genome::DataSource::Transcripts {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton'],
};

sub delimiter {
    return ",";
}

sub column_order {
    return [qw(
        transcript_id
        gene_id
        transcript_start
        transcript_stop
        transcript_name
        transcript_status
        strand
        chrom_name
        species 
        source 
        version
        gene_name
        transcript_error
        coding_region_start
        coding_region_stop
        amino_acid_length
    )]
}

sub constant_values { [qw/ data_directory reference_build_id /] }; 

sub sort_order {
    return [qw(chrom_name transcript_start transcript_stop transcript_id)];
}

sub skip_first_line {
    return 0;
}

sub required_for_get { ['data_directory', 'reference_build_id'] }

sub file_resolver {
    my ($data_directory, $reference_build_id) = @_;

    my $path = "$data_directory/transcripts.csv";

    return $path;
}

1;

