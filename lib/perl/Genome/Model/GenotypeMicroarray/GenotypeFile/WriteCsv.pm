package Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Writer;
use Genome::Utility::IO::SeparatedValueWriter;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv { 
    is => 'Genome::Utility::IO::SeparatedValueWriter',
};

sub create {
    my ($class, %params) = @_;

    delete $params{sample_name};
    $params{separator} = "\t" if not $params{separator} or $params{separator} =~ /^tab$/i;
    $params{print_headers} = delete $params{headers} if exists $params{headers};
    $params{in_place_of_null_value} = 'NA';
    $params{ignore_extra_columns} = 1;

    my $fields = delete $params{fields};
    if ( $fields ) {
        $params{headers} = [ split(',', $fields) ];
    }
    else {
        $params{headers} = [qw/ chromosome position alleles reference id sample_name log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /];
    }

    return $class->SUPER::create(%params);
}

sub write {
    my ($self, $entry) = @_;

    $self->write_one($entry);
    return;

    my %genotype;
    #$self->write_one(\%genotype);
    
    return 1;
}

1;

