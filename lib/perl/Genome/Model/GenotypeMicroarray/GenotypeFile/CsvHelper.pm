package Genome::Model::GenotypeMicroarray::GenotypeFile::CsvHelper;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::CsvHelper { 
    is => 'UR::Singleton',
};

sub column_names {
    return [qw/ chromosome position alleles reference id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /];
}

1;

