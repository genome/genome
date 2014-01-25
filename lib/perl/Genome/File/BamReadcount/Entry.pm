package Genome::File::BamReadcount::Entry;

use Genome;
use Genome::File::BamReadcount::LibraryMetrics;
use Genome::File::BamReadcount::AlleleMetrics;

use strict;
use warnings;

sub new {
    my ($class, $string) = @_;
    my $self;
    my ($chromosome, $position, $ref_base, $depth, $rest) = 
        $string =~ /^
            (.+)\t
            (\d+)\t
            ([AaCcGgTtRrYyKkMmSsWwBbDdHhVvNn])\t
            (\d+)\t
            (.+)$/x;
    unless(defined $chromosome and defined $position and defined $ref_base and defined $depth and defined $rest) {
        die "Error parsing out Entry\n";
    }

    $self = {
        _chromosome => $chromosome,
        _position => $position,
        _ref_base => $ref_base,
        _depth => $depth,
        _source_string => $string,
        _allele_metrics => undef,
        _alleles => undef,
        _allele_map => undef,
        _libraries => undef,
    };

    my @library_strings = ( $string =~ /(\S+\t\{.+?\})/g );
    my @libraries = map { Genome::File::BamReadcount::LibraryMetrics->new($_) or die "Unable to create LibraryMetrics object from $_" } @library_strings;

    unless(@libraries) {
        #must not have been run in per_lib mode
        my @allele_metrics = map { Genome::File::BamReadcount::AlleleMetrics->new($_) or die "Unable to create AlleleMetrics object from $_"; } split /\t/, $rest;
        my @alleles = map { $_->allele } @allele_metrics;
        my %allele_map;
        @allele_map{@alleles} = @allele_metrics;
        $self->{_allele_metrics} =  \%allele_map;
        $self->{_alleles} = \@alleles;
    }
    else {
        $self->{_libraries} = \@libraries;
    }
    bless $self, $class;
}

sub chromosome {
    my ($self) = @_;
    return $self->{_chromosome};
}

sub position {
    my ($self) = @_;
    return $self->{_position};
}

sub ref_base {
    my ($self) = @_;
    return $self->{_ref_base};
}

sub depth {
    my ($self) = @_;
    return $self->{_depth};
}

sub libraries {
    my ($self) = @_;
    if(defined $self->{_libraries}) {
        return @{$self->{_libraries}};
    }
    else {
        return undef;
    }
}

sub num_libraries {
    my ($self) = @_;
    if(defined $self->libraries) {
        return scalar($self->libraries);
    }
    else {
        return 0;
    }
}

sub has_per_library {
    my ($self) = @_;
    return defined $self->{_libraries};
}

sub metrics_for {
    my ($self, $allele) = @_;
    if($self->has_per_library) {
        #TODO either add whole entry metrics to the executable or try to form some sort of perl summary here
        die "Please grab metrics through the library objects instead of the Entry directly.\n";
    }
    else {
        if(exists $self->{_allele_metrics}{$allele}) {
            return $self->{_allele_metrics}{$allele};
        }
        else {
            die "Can't retrieve metrics for allele '$allele'\n";
        }
    }
}

sub alleles {
    my ($self) = @_;
    if($self->has_per_library) {
        #TODO either add whole entry support to the executable or try to form some sort of perl summary here
        die "Please grab alleles through the library objects instead of the Entry directly.\n";
    }
    else {
        return @{$self->{_alleles}};
    }
}


1;
