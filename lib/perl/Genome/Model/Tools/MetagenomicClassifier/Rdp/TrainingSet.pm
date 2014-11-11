package Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet;

use strict;
use warnings;

use Genome;

require Bio::Taxon;
require List::MoreUtils;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet {
    id_by => {
        training_path => {
            is => 'Text',
            doc => 'Training path to inspect.'
        },
    },
};

sub valid_training_sets {
    return (qw/ 4 6 9 10 broad /);
}

sub base_training_path {
    return "/gsc/scripts/share/rdp/";
}

sub create_from_training_set_name {
    my ($class, $training_set_name) = @_;

    die $class->error_message('No training set given to get training path!') if not defined $training_set_name;
    if ( List::MoreUtils::none { $training_set_name eq $_ } $class->valid_training_sets ) {
        die $class->error_message("Invalid training set ($training_set_name) given to create_from_training_set_name!");
    }

    my $training_path = $class->base_training_path.'/'.$training_set_name;
    if ( not -d $training_path ) {
        die $class->error_message("Training path does not exist: $training_path");
    }

    return $class->create(training_path => $training_path);
}

sub classifier_properties_path {
    my $self = shift;

    my $classifier_properties_path = $self->training_path.'/rRNAClassifier.properties';
    if ( not -s $classifier_properties_path ) {
        die $self->error_message('No rRNAClassifier.properties in training path! '.$self->training_path);
    }

    return $classifier_properties_path;
}


sub get_genera {
    my $self = shift;
    unless ($self->{_genera}) {
        $self->_build_taxonomy();
    }
    return $self->{_genera};
}

sub get_taxonomy {
    my $self = shift;
    unless ($self->{_taxonomy}) {
        $self->_build_taxonomy();
    }
    return $self->{_taxonomy};
}

sub _get_taxonomy_path {
    my $self = shift;
    my $path = $self->training_path;
    my $taxonomy_path = $path."/bergeyTrainingTree.xml"; 
    unless (-e $taxonomy_path) {
        die __PACKAGE__.": can't find bergeyTrainingTree.xml in $path";
    }
    return $taxonomy_path;
}

sub _build_taxonomy {
    my $self = shift;

    my $in = Genome::Sys->open_file_for_reading($self->_get_taxonomy_path);

    my %taxons;
    my @genera;

    while (my $line = <$in>) {
        if ($line =~ /TreeNode/) {
            chomp $line;
            $line =~ s/TreeNode[ >]//g;
            $line =~ s/[<>\/]//g;
            $line =~ s/'//g;
            $line =~ s/" /:/g;
            $line =~ s/"//g;
            my %attributes = split /[:=]/,$line;
            my $taxon =   new Bio::Taxon(
                             -id        => $attributes{taxid},
                             -rank      => $attributes{rank},
                         );
            $taxon->node_name($attributes{name});
            $taxons{$taxon->id} = $taxon;
            my $parent = $taxons{$attributes{parentTaxid}};
            if ($parent && $attributes{taxid} != $attributes{parentTaxid}) {
                $parent->add_Descendent($taxon);
            }
            if ($taxon->rank eq 'genus') {
                push @genera, $taxon;
            }
        }
    }

    $self->{_taxonomy} = $taxons{0};
    $self->{_genera} = \@genera;
}

1;

