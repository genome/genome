package Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet;

use strict;
use warnings;

use Genome;

require Bio::Taxon;
require File::Spec;
require List::MoreUtils;

class Genome::Model::Tools::MetagenomicClassifier::Rdp::TrainingSet {
    id_by => {
        path => {
            is => 'Text',
            doc => 'Training path to inspect.'
        },
    },
};

sub valid_set_names {
    return (qw/ 4 6 9 10 broad /);
}

sub base_path {
    return File::Spec->join(File::Spec->rootdir, 'gsc', 'scripts', 'share', 'rdp'); # FIXME should be ENV
}

sub path_for_set_name {
    my ($class, $set_name) = @_;

    die $class->error_message('No set name given to get path!') if not defined $set_name;
    if ( List::MoreUtils::none { $set_name eq $_ } $class->valid_set_names ) {
        die $class->error_message("Invalid set name! $set_name");
    }

    return File::Spec->join(base_path(), $set_name);
}

sub classifier_properties_path {
    my $self = shift;

    my $classifier_properties_path = File::Spec->join($self->path, 'rRNAClassifier.properties');
    if ( not -s $classifier_properties_path ) {
        die $self->error_message('No classifier properties (rRNAClassifier.properties) in training path! '.$self->path);
    }

    return $classifier_properties_path;
}

sub taxonomy_path {
    my $self = shift;
    my $path = $self->path;
    my $taxonomy_path = File::Spec->join($path, "bergeyTrainingTree.xml"); 
    unless (-e $taxonomy_path) {
        die $self->error_message("No taxonomy XML (bergeyTrainingTree.xml) in training path! $path");
    }
    return $taxonomy_path;
}


sub get_genera {
    my $self = shift;
    unless ($self->{_genera}) {
        my $build_ok = $self->_build_taxonomy();
        return if not $build_ok;
    }
    return $self->{_genera};
}

sub get_taxonomy {
    my $self = shift;
    unless ($self->{_taxonomy}) {
        my $build_ok = $self->_build_taxonomy();
        return if not $build_ok;
    }
    return $self->{_taxonomy};
}

sub _build_taxonomy {
    my $self = shift;

    my $in = Genome::Sys->open_file_for_reading($self->taxonomy_path);

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
            $attributes{name} =~ s/&quot;//g;
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

