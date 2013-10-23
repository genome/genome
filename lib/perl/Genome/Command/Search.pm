package Genome::Command::Search;

use strict;
use warnings;


class Genome::Command::Search {
    is => 'Command::V2',

    has => [
        target => {
            is => 'Text',
            shell_args_position => 1,
            is_many => 1,
            doc => 'The string, id, or partial id to search for.',
        },
    ],
};


sub execute {
    my $self = shift;

    $self->show($self->get_docs);

    return 1;
}


sub show {
    my ($self, $docs) = @_;
    my $type = $docs->[0]->{type};

    my @objects = $self->get_objects($docs);
    $self->show_objects($type, @objects);
}


my %OBJECT_GETTERS = (
    'build' => 'Genome::Model::Build',
#    "illumina_run"
#    "imported_instrument_data"
#    "individual"
#    "instrument data"
#    "library"
#    "mail"
    "model - ClinSeq" => 'Genome::Model',
    "model - alignment" => 'Genome::Model',
    "model - convergence" => 'Genome::Model',
    "model - lane_qc" => 'Genome::Model',
    "model - microarray" => 'Genome::Model',
    "model - other" => 'Genome::Model',
    "model - rna" => 'Genome::Model',
    "model - somatic" => 'Genome::Model',
    "modelgroup" => 'Genome::ModelGroup',
#    "population_group"
#    "processing_profile"
#    "project"
#    "sample"
#    "solexa_instrument_data"
#    "taxon"
#    "user"
#    "wiki-page"
#    "work-order"
);

sub get_objects {
    my ($self, $docs) = @_;

    my $end = $self->find_different_getter($docs);

    my @useful_docs = (@$docs)[0..$end];
    my @ids = map {$_->{object_id}} @useful_docs;

    my $getter = $OBJECT_GETTERS{$docs->[0]->{type}};
    return $getter->get(id => \@ids);
}

sub find_different_getter {
    my ($self, $docs) = @_;
    my $len = scalar @$docs;

    my $getter = $OBJECT_GETTERS{$docs->[0]->{type}};
    for (my $i = 0; $i < $len; $i++) {
        if (exists $OBJECT_GETTERS{$docs->[$i]->{type}}) {
            my $this_getter = $OBJECT_GETTERS{$docs->[$i]->{type}};
            if ($getter ne $this_getter) {
                return $i - 1;
            }
        } else {
            return $i - 1;
        }
    }
    return $len - 1;
}

my %OBJECT_SHOWERS = (
    'build' => 'Genome::Command::Search::Build',
#    "illumina_run"
#    "imported_instrument_data"
#    "individual"
#    "instrument data"
#    "library"
#    "mail"
    "model - ClinSeq" => 'Genome::Command::Search::Model',
    "model - alignment" => 'Genome::Command::Search::Model',
    "model - convergence" => 'Genome::Command::Search::Model',
    "model - lane_qc" => 'Genome::Command::Search::Model',
    "model - microarray" => 'Genome::Command::Search::Model',
    "model - other" => 'Genome::Command::Search::Model',
    "model - rna" => 'Genome::Command::Search::Model',
    "model - somatic" => 'Genome::Command::Search::Model',
    "modelgroup" => 'Genome::Command::Search::ModelGroup',
#    "population_group"
#    "processing_profile"
#    "project"
#    "sample"
#    "solexa_instrument_data"
#    "taxon"
#    "user"
#    "wiki-page"
#    "work-order"
);

sub show_objects {
    my ($self, $type, @objects) = @_;

    my $shower = $OBJECT_SHOWERS{$type};

    if (scalar(@objects) > 1) {
        $shower->display_many(\@objects);
    } else {
        $shower->display_single($objects[0]);
    }
}

sub get_docs {
    my $self = shift;

    my $content = $self->get_content;
    return $content->{response}{docs};
}

sub get_content {
    my $self = shift;

    my $response = Genome::Search->search(join(' ', $self->target));
    unless (defined($response)) {
        die "Invalid response from Search";
    }

    $self->validate_response_content($response->content);

    return $response->content;
}


my $SCORE_WARNING_THRESHOLD = 5;
sub validate_response_content {
    my ($self, $content) = @_;


    if ($content->{response}{numFound} <= 0) {
        die "No results found";
    }

    if ($content->{response}{maxScore} < $SCORE_WARNING_THRESHOLD) {
        warn sprintf("Warning: Low search score (%s < %s).",
            $content->{response}{maxScore}, $SCORE_WARNING_THRESHOLD);
    }
}


1;
