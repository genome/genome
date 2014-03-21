package Genome::Model::Tools::Graph::MutationDiagram::DomainProvider;

use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Graph::MutationDiagram::DomainProvider {
    has => [
        custom_domains => {
            is => "String",
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);

    if ($self->custom_domains) {
        my @domain_specification = split(',', $self->custom_domains);

        my @custom_domains;
        while(@domain_specification) {
            my %domain = (type => "CUSTOM");
            @domain{qw(name start end)} = splice @domain_specification,0,3;
            push @custom_domains, \%domain;
        }

        $self->{_custom_domains} = \@custom_domains;
    }
    return $self;
}

sub get_domains {
    my $self = shift;
    if ($self->{_custom_domains}->[0]) {
        return @{$self->{_custom_domains}};
    }
    return;
}

sub get_amino_acid_length {
    my $self = shift;
    die $self->error_message("Must implement get_amino_acid_length in child class");
}

1;

