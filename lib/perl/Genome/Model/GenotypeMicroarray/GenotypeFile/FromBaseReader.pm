package Genome::Model::GenotypeMicroarray::GenotypeFile::FromBaseReader;

use strict;
use warnings;

use Data::Dumper 'Dumper';

our %source_types_and_delimiters = (
    build => "\t",
    instrument_data => ',',
);

sub create {
    my ($class, %params) = @_;

    my $self = bless(\%params, $class);

    my $source_type = $self->source_type;
    if ( not $self->{$source_type} ) {
        print STDERR "No $source_type given to open genotype file reader!";
        return;
    }

    $self->{delimiter} = $source_types_and_delimiters{$source_type};

    my $open_genotype_fh_ok = $self->_open_genotype_fh;
    return if not $open_genotype_fh_ok;

    my $resolve_headers = $self->_resolve_headers;
    return if not $resolve_headers;

    return $self;
}

sub _open_genotype_fh {
    my $self = shift;

    my $genotype_file = $self->get_genotype_file;
    if ( not $genotype_file or not -s $genotype_file ) {
        print STDERR 'No original genotype file for build! '.Dumper($self->{ $self->source_type });
        return;
    }

    my $fh = eval { Genome::Sys->open_file_for_reading($genotype_file) };
    if (!$fh or $@) {
        $self->error_message("Can't open file $genotype_file for reading: $@");
        return;
    }

    $self->{genotype_fh} = $fh;

    return 1;
}

sub read {
    my $self = shift;

    my $line = $self->{genotype_fh}->getline;
    return if not $line;
    chomp $line;

    my %genotype;
    @genotype{ @{$self->{headers}} } = split($self->{delimiter}, $line);

    return \%genotype;
}

1;

