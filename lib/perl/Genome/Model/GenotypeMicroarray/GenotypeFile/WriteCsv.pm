package Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Writer;
use Genome::Utility::IO::SeparatedValueWriter;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriteCsv { 
    is => 'Genome::Utility::IO::SeparatedValueWriter',
    has => {
        _format_types => { is => 'Array', },
        _sample_name => { is => 'Text', },
    },
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

    my $header = delete $params{header};
    if ( not $header ) {
        $class->error_message('No header given!');
        return;
    }

    my $self = $class->SUPER::create(%params);
    return if not $self;

    $self->{_sample_name} = ($header->sample_names)[0];
    $self->{_format_types} = [ keys %{$header->format_types} ] if $header->format_types;

    return $self;
}

sub write {
    my ($self, $entry) = @_;

    my %genotype = (
        id => $entry->{identifiers}->[0],
        chromosome => $entry->{chrom},
        position => $entry->{position},
        reference => $entry->{reference_allele},
        sample_name => $self->_sample_name,
    );

    for my $format_type ( @{$self->_format_types} ) {
        my $value = $entry->sample_field(0, $format_type);
        $value = 'NA' if not defined $value or $value eq '.';
        $genotype{ Genome::Model::GenotypeMicroarray->format_name_for_id($format_type) } = $value;
    }
    @genotype{qw/ allele1 allele2 /} = split(//, $genotype{alleles});

    $self->write_one(\%genotype);

    return 1;
}

1;

