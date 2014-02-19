package Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv;

use strict;
use warnings;

use Genome;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv { 
    is => 'Genome::Utility::IO::SeparatedValueWriter',
    has => {
        ignore_extra_columns => { value => 1, },
    },
};

sub create {
    my ($class, %params) = @_;

    my $process_params = $class->_process_params(\%params);
    return if not $process_params;

    my $self = $class->SUPER::create(%params);
    return if not $self;

    return $self;
}

sub _process_params {
    my ($class, $output_params) = @_;

    $output_params->{separator} = "\t" if not $output_params->{separator} or $output_params->{separator} =~ /^tab$/i;
    $output_params->{print_headers} = delete $output_params->{headers} if exists $output_params->{headers};

    my $fields = delete $output_params->{fields};
    if ( $fields ) {
        $output_params->{headers} = [ split(',', $fields) ];
    }
    else {
        $output_params->{headers} = [qw/ chromosome position alleles reference id sample_id log_r_ratio gc_score cnv_value cnv_confidence allele1 allele2 /];
    }

    return 1;
}

1;

