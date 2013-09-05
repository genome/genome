package Genome::DruggableGene::GeneNameGroup;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::GeneNameGroup {
    is => 'UR::Object',
    table_name => 'dgidb.gene_name_group',
    schema_name => 'dgidb',
    data_source => 'Genome::DataSource::Dgidb',

    id_generator => '-uuid',
    id_by => [
        id => {is => 'Text'},
    ],
    has => [
        display_name => {
            is => 'Text',
            calculate_from => ['name'],
            calculate => q| $name |,
        },
        name => { is => 'Text' },
        bridges => {
            is => 'Genome::DruggableGene::GeneNameGroupBridge',
            reverse_as => 'group',
            is_many => 1,
        },
        genes  => {
            is => 'Genome::DruggableGene::GeneNameReport',
            via => 'bridges',
            to => 'gene',
            is_many => 1,
        },
    ],
    doc => 'Group of likely synonymous genes',
};

sub consume {
    my $self = shift;
    my @groups = @_;
    for my $group (@groups){
        for my $bridge($group->bridges){
             $self->add_bridge(gene_id => $bridge->gene_id);
        }
        $group->delete;
    }
}

sub delete {
    my $self = shift;
    for($self->bridges) {
        $_->delete;
    }
    return $self->SUPER::delete();
}

1;
