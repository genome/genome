package Genome::Model::Tools::MetagenomicClassifier::Rdp::ListGenera;

use strict;
use warnings;

require Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::ListGenera {
    is => 'Command',
    has_optional => {
        training_set_name => {
            type => 'String',
            valid_values => [ Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->valid_set_names ],
            doc => 'Name of training set to inspect. Uses the base training set path: '
                    . Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->base_path,
        },
        training_set_path => {
            type => 'String',
            doc => 'Path to training set.',
        },
    },
    has_optional_transient => {
        training_set => { is => 'Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet', },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( not defined $self->training_set and not defined $self->training_set_path ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ path set_name /],
                desc => 'Neither training set name nor path is specified! Please indicate which option to use.', 
            )
        );
    }

    if ( not -d $self->training_set_path ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ path /],
                desc => 'Traning set path does not exist! '.$self->training_set_path, 
            )
        );
    }

    return;
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    if ( defined $self->training_set_name ) {
        my $training_set_path = Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->path_for_set_name(
            $self->training_set_name
        );
        $self->training_set_path($training_set_path);
    }

    return $self;
}

sub execute {
    my $self = shift;

    my $training_set = Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->create(
        path => $self->training_set_path,
    );
    return if not $training_set;
    $self->training_set($training_set);

    foreach my $genus ( @{$training_set->get_genera} ) {
        print _to_string($genus) . "\n";
    }

    return 1;
}

sub _to_string {
    my $genus = shift;
    my @ancestors;

    my @ranks = ('no rank', 'domain',  'phylum',  'class', 'order', 'family', 'genus');
    my $current_taxon = $genus;
    while (defined $current_taxon) {
        unshift @ancestors, $current_taxon;
        $current_taxon = $current_taxon->ancestor
    }

    my $taxonomy_string = "";
    my $first = 1;
    my $rank = 1;

    my $root = shift @ancestors;
    $taxonomy_string = $root->node_name.";";
    foreach my $taxon (@ancestors) {
        next if lc($taxon->rank) =~ /sub|super/;
        if ($first) {
            $first = 0;
        }
        else {
            $taxonomy_string.= ";";
        }

        while (lc($taxon->rank) !~ $ranks[$rank]) {
            $taxonomy_string .= "unknown_".$ranks[$rank++].";"; 
        }
        $taxonomy_string .= $taxon->node_name;
        $rank++;
    }

    return $taxonomy_string;
}

#< HELP >#
sub help_brief {
    "rdp training-set taxonomy lister",
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
genome-model tools metagenomic-classifier rdp list-genera   
EOS
}

1;

