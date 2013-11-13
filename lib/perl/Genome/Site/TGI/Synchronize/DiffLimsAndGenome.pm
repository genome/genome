package Genome::Site::TGI::Synchronize::DiffLimsAndGenome;

use strict;
use warnings;

use Genome;

use Genome::Site::TGI::Synchronize::Classes::Dictionary;
use Set::Scalar;

class Genome::Site::TGI::Synchronize::DiffLimsAndGenome {
    is => 'Command::V2',
    has => [
        entity_name => {
            is => 'Text',
            valid_values => [ Genome::Site::TGI::Synchronize::Classes::Dictionary->entity_names ],
            doc => 'Entity to diff between LIMS and Genome.',
        },
    ],
    has_optional => [
        print_diffs => { is => 'Boolean', default_value => 1, doc => 'Print diffs of IDS to STDOUT.', },
    ],
    has_optional_transient => [
        in_lims_not_genome => { is => 'Set::Scalar', },
        in_genome_not_lims => { is => 'Set::Scalar', },
    ],
    doc => 'Generate diffs of IDs in LIMS and Genome for an entity'
};

sub lims_class {
    return Genome::Site::TGI::Synchronize::Classes::Dictionary->lims_class_for_entity_name($_[0]->entity_name);
}

sub genome_class {
    return Genome::Site::TGI::Synchronize::Classes::Dictionary->genome_class_for_entity_name($_[0]->entity_name);
}

sub execute {
    my $self = shift;

    my $entity_name = $self->entity_name;
    $self->status_message("Diff LIMS and Genome $entity_name...");

    my $lims_class = $self->lims_class;
    my $lims_ids = $self->_create_id_set_for_class($lims_class);

    my $genome_class =$self->genome_class;
    my $genome_ids = $self->_create_id_set_for_class($genome_class);

    my $in_lims_not_genome = $lims_ids->difference($genome_ids);
    $self->in_lims_not_genome($in_lims_not_genome);

    my $in_genome_not_lims = $genome_ids->difference($lims_ids);
    $self->in_genome_not_lims($in_genome_not_lims);

    $self->_print_diffs if $self->print_diffs;

    $self->status_message("Diff LIMS and Genome $entity_name...done");
    return 1;
}

sub _create_id_set_for_class {
    my ($self, $class) = @_;
    $self->status_message("Getting IDs for $class...");

    my $iterator = $class->create_iterator;
    my $set = Set::Scalar->new();
    while ( my $obj = $iterator->next ) {
        $set->insert($obj->id);
    };
    $self->status_message("Found ".scalar(@$set)." $class IDs.");

    return $set;
}

sub _print_diffs {
    my $self = shift;

    my $entity_name = $self->entity_name;
    print join(
        "\n",
        uc($entity_name).'S IN LIMS NOT GENOME:',
        $self->in_lims_not_genome->elements,
        uc($entity_name).'S IN GENOME NOT LIMS:',
        $self->in_genome_not_lims->elements,
    )."\n";

    return 1;
}

1;

