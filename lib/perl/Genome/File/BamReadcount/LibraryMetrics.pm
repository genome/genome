package Genome::File::BamReadcount::LibraryMetrics;

use Genome::File::BamReadcount::AlleleMetrics;
use Genome;

use strict;
use warnings;

sub new {
    my ($class, $string) = @_;
    my $self;
    my ($name, $allele_metric_string) = $string =~ /^
        (\S+?) #library name. Can't have any whitespace
        \t\{\t #tab bracket tab sets off the allele metrics
        (.+?)  #allele metrics
        \t\}$/x;#ends with closing bracket
    unless(defined $name and defined $allele_metric_string) {
        die "Unable to create LibraryMetrics from string '$string'\n";
    }
    my @allele_metrics = map { Genome::File::BamReadcount::AlleleMetrics->new($_) or die "Unable to create AlleleMetrics object from $_"; } split /\t/, $allele_metric_string;
    $self = { _name => $name };
    my @alleles = map { $_->allele } @allele_metrics;
    my %allele_map;
    @allele_map{@alleles} = @allele_metrics;
    $self->{_alleles} = \@alleles;
    $self->{_allele_metrics} = \%allele_map;

    bless $self, $class;
    return $self;
}

sub alleles {
    my ($self) = @_;
    return @{$self->{_alleles}};
}

sub depth {
    my ($self) = @_;
    my $total_count = 0;
    for my $allele ($self->alleles) {
        $total_count += $self->{_allele_metrics}{$allele}->count;
    }
    return $total_count;
}

sub metrics_for {
    my ($self, $allele) = @_;
    if(exists $self->{_allele_metrics}{$allele}) {
        return $self->{_allele_metrics}{$allele};
    }
    else {
        die "Allele '$allele' not found in LibraryEntry.";
    }
}

sub name {
    my ($self) = @_;
    return $self->{_name};
}

1;
