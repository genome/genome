package Genome::Utility::MetagenomicClassifier::Rdp::TrainingSet;

use strict;
use warnings;

use Data::Dumper;
require Bio::Taxon;
require Genome::Sys;
require IO::File;

sub create {
    my ($class, %params) = @_;
    
    unless (exists $params{path}) {
        die __PACKAGE__." create requires path";
    }
    my $self = bless \%params, $class;

    return $self;
}

sub get_genera {
    my $self = shift;
    unless ($self->{genera}) {
        $self->_build_taxonomy();
    }
    return $self->{genera};
}

sub get_taxonomy {
    my $self = shift;
    unless ($self->{taxonomy}) {
        $self->_build_taxonomy();
    }
    return $self->{taxonomy};
}

sub _get_taxonomy_path {
    my $self = shift;
    my $path = $self->path;
    my $taxonomy_path = $path."/bergeyTrainingTree.xml"; 
    unless (-e $taxonomy_path) {
        die __PACKAGE__.": can't find bergeyTrainingTree.xml in $path";
    }
    return $taxonomy_path;
}

sub path {
    my $self = shift;
    return $self->{path};
}

sub _build_taxonomy {
    my $self = shift;
    my $in = new IO::File($self->_get_taxonomy_path);

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

    $self->{taxonommy} = $taxons{0};
    $self->{genera} = \@genera;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Utility/MetagenomicClassifier/Rdp/TrainingSet.pm $
#$Id: $
