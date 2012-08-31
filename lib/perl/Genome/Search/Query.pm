
package Genome::Search::Query;

use strict;
use warnings;

our $PAGE_SIZE = 50;

class Genome::Search::Query {
    is => 'UR::Value',
    id_by => [
        query => { is => 'Text' },
        page => { is => 'Number' },
        fq   => { is => 'Text' }
    ],
    has => [
        results => {
            is => 'Genome::Search::Result',
            reverse_as => 'query',
            is_many => 1,
        },
        result_subjects => {
            is_many => 1,
            via => 'results',
            to => 'subject'
        },
        page_size => {
            is_constant => 1,
            is_class_wide => 1,
            value => $PAGE_SIZE
        },
        total_found => {
            is => 'Number'
        },
        page_count => {
            is => 'Number'
        }
    ]
};

sub get {
    my $class = shift;
    if (@_ > 1 && @_ % 2 == 0) {
        my %args = (@_);
        $args{'fq'} ||= '';
        $args{'page'} = 1 unless exists $args{page};
        @_ = %args;
    }

    return $class->SUPER::get(@_);
}

sub total_found {
    my $self = shift;

    unless (exists $self->{executed} && $self->{executed}) {
        $self->execute;
    }
    $self->__total_found;
}

sub results {
    my $self = shift;

    unless (exists $self->{executed} && $self->{executed}) {
        $self->execute;
    }

    return $self->__results;
}

sub result_objects {
    my $self = shift;

    unless (exists $self->{executed} && $self->{executed}) {
        $self->execute;
    }

    return $self->__result_objects;
}

sub execute {
    my $self = shift;
    $self->{executed} = 1;

    my $response = Genome::Search->search(
        $self->query,
        {
            rows  => $PAGE_SIZE,
            start => $PAGE_SIZE * ( $self->page - 1 )
        }
    );

    foreach my $doc ($response->docs) {
        my $subject_class = $doc->value_for('class');
        my $subject_id = $doc->value_for('object_id');

        Genome::Search::Result->create(
            query_string => $self->query,
            page => $self->page,
            subject_class_name => $subject_class,
            subject_id => $subject_id,
            fq => $self->fq(),
        );
    }

    my $pages = $response->content->{response}->{numFound} / $self->page_size;

    if ($pages == int($pages)) {
        $self->page_count($pages);
    } else {
        $self->page_count(1 + int($pages));
    }
   
    $self->__total_found($response->content->{response}->{numFound});
}

