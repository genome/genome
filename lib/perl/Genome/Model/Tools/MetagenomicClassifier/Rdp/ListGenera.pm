package Genome::Model::Tools::MetagenomicClassifier::Rdp::ListGenera;

use strict;
use warnings;

require Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::ListGenera {
    is => 'Command',
    has_optional => {
        training_set => {
            type => 'String',
            valid_values => [ Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->valid_training_sets ],
            doc => 'Name of training set (broad)',
        },
        training_path => {
            type => 'String',
            doc => 'Path to training set.',
        },
    },
    has_optional_transient => {
        _training_set => { is => 'Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet', },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    if ( $self->training_set and $self->training_path ) {
        return (
            UR::Object::Tag->create(
                type => 'invalid',
                properties => [qw/ training_path training_set /],
                desc => 'Both training_set and training_path are specified! Please only include one option.', 
            )
        );
    }

    my $training_set;
    if ( $self->training_set ) {
        $training_set = Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->create_from_training_set_name(
            $self->training_set
        );
        $self->training_path( $training_set->training_path );
    }
    else {
        $training_set = Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet->create(
            training_path => $self->training_path,
        );
    }
    $self->_training_set($training_set);

    return;
}

sub execute {
    my $self = shift;

    foreach my $genus ( @{$self->_training_set->get_genera} ) {
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

