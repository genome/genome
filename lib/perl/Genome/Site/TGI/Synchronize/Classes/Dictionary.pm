package Genome::Site::TGI::Synchronize::Classes::Dictionary;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Site::TGI::Synchronize::Classes::Dictionary { 
    is => 'UR::Singleton',
};

our @lims_classes = (
    map { 
        'Genome::Site::TGI::Synchronize::Classes::'.$_
    } (qw/ 
        OrganismTaxon OrganismIndividual PopulationGroup OrganismSample LibrarySummary
        LimsProject LimsProjectSample LimsProjectInstrumentData
        RegionIndex454 IndexIllumina Genotyping
    /)
);

sub entity_names {
    return map { $_->entity_name } grep { $_->can('entity_name') } @lims_classes;
}

sub lims_class_for_entity_name {
    my ($self, $entity_name) = @_;

    Carp::confess( $self->error_message('No entity name given to get LIMS class!') )if not $entity_name;

    for my $lims_class ( @lims_classes ) {
        next if not $lims_class->can('entity_name');
        return $lims_class if $lims_class->entity_name eq $entity_name;
    }

    $self->error_message("No LIMS class for entity! $entity_name");
    return;
}

sub genome_class_for_entity_name {
    my ($self, $entity_name) = @_;

    my $lims_class = $self->lims_class_for_entity_name($entity_name);
    return if not $lims_class;

    return $lims_class->genome_class_for_create;
}

1;

