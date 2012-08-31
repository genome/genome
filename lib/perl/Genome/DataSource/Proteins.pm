package Genome::DataSource::Proteins;

use strict;
use warnings;

use Genome;

class Genome::DataSource::Proteins {
    is => [ 'UR::DataSource::FileMux', 'UR::Singleton'],
};

sub delimiter {
    return ",";
}

sub column_order {
    return [qw( 
    protein_id 
    transcript_id 
    protein_name 
    amino_acid_seq 
    species 
    source 
    version
    )];
}

sub sort_order {
    return [qw( transcript_id protein_id )];
}

sub skip_first_line {
    return 0;
}

sub constant_values { [qw/ data_directory reference_build_id /] };
sub required_for_get { [qw( transcript_id data_directory reference_build_id)] }

sub file_resolver {
    my($composite_id, $data_directory) = @_;

    my $meta = Genome::Transcript->__meta__;

    my ($chrom, $start, $stop, $species, $source, $version, $transcript_id) = $meta->resolve_ordered_values_from_composite_id($composite_id);
    
    my $thousand = int($transcript_id / 1000);
    my $path = "$data_directory/proteins/proteins_" . $thousand . ".csv";
    return $path;
}

1;

