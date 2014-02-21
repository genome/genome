package Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf;

use strict;
use warnings;

use Genome;

use Genome::File::Vcf::Writer;

class Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf { 
    is => 'UR::Object',
    has => {
        output => { is => 'Text', },
        header => { is => 'Genome::File::Vcf::Header', },
    },
    has_transient_optional => {
        _writer => { is => 'Genome::File::Vcf::Writer', },
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

    $self->_writer(
        Genome::File::Vcf::Writer->new($self->output, $self->header)
    );

    return $self;
}

sub write_one {
    my ($self, $genotype) = @_;
    $self->_writer->write($genotype);
    return 1;
}

sub close {
    return $_[0]->_writer->close;
}

1;

