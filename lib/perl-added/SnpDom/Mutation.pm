package SnpDom::Mutation;

use strict;
use warnings;
use Carp;

sub new {
    my ($class,%args) = @_;
    my $self = { "_gene" => undef,
                 "_transcript" => undef,
                 "_gt" => undef,
                 "_mutations" => [ ],
                 "_domains" => undef,
                };
    if(exists($args{gene})) {
        $self->{_gene} = $args{gene};
    }
    if(exists($args{'gt'})) {
        $self->{_gt} = $args{'gt'};
    }
    bless($self,$class);
    return $self;
}

sub add_mutation {
    my ($self,$mutation_string) = @_;
    push(@{$self->{_mutations}},$mutation_string);
    return 1;
}

sub add_domain {
    my ($self,$aachange,$domain) = @_;
    push(@{$self->{_domains}->{$aachange}},$domain);
    return 1;
}

sub get_domain {
    my ($self,$aachange) = @_;
    if(exists($self->{_domains}->{$aachange}) ) {
        return $self->{_domains}->{$aachange};
    }
    else {
        return;
    }
    return;
}

sub get_all_mutations {
    my ($self) = @_;
    #return @{$self->{_mutations}};
    return $self->{_mutations};
}

sub get_count_of_mutations {
    my ($self) = @_;
    return $#{$self->{_mutations}} + 1;
}

sub get_name {
    my ($self) = @_;
    return $self->{_gt};
}

sub set_name {
    my ($self,%args) = @_;
    if(exists($args{transcript}))  {
        $self->{_transcript} = $args{transcript};
    }
    if(exists($args{'gt'}))  {
        $self->{_gt} = $args{'gt'};
    }
    if(exists($args{gene}))  {
        $self->{_gene} = $args{gene};
    }

    return 1;
}

1;
