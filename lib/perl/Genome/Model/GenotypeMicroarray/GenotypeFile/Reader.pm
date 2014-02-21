package Genome::Model::GenotypeMicroarray::GenotypeFile::Reader;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Entry;

class Genome::Model::GenotypeMicroarray::GenotypeFile::Reader { 
    is => 'UR::Object',
    has => {
        reader => { is => 'UR::Object', },
    },
    has_optional_transient => {
        _read_from_reader => { is => 'Code', },
    },
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my @errors = $self->__errors__;
    if ( @errors ) {
        Carp::confess( join ("\n", map { $_->__display_name__ } @errors) );
    }

    $self->_resolve_read_sub;

    return $self;
}

sub _resolve_read_sub {
    my $self = shift;

    my $reader_class = $self->reader->class;
    if ( $reader_class eq 'Genome::File::Vcf::Reader' ) {
        $self->_read_from_reader(
            sub {
                my $self = shift;
                return $self->reader->next;
            }
        );
    }
    elsif ( $reader_class =~ /^Genome::Model::GenotypeMicroarray::GenotypeFile::ReadTsv/ ) {
        $self->_read_from_reader(
            sub {
                my $self = shift;

                my $genotype = $self->next;
                return if not $genotype;

                $genotype->{alleles} = $genotype->{allele1}.$genotype->{allele2};

                return $genotype;
            }
        );
    }
    else {
        Carp::confess('Unknown reader! '.$self->reader);
    }

    return $self->_read_from_reader;
}

1;

