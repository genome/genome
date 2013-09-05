package Genome::DruggableGene::DrugNameGroup;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::DrugNameGroup {
    is => 'UR::Object',
    table_name => 'dgidb.drug_name_group',
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
            is => 'Genome::DruggableGene::DrugNameGroupBridge',
            reverse_as => 'group',
            is_many => 1,
        },
        drugs  => {
            is => 'Genome::DruggableGene::DrugNameReport',
            via => 'bridges',
            to => 'drug',
            is_many => 1,
        },
    ],
    doc => 'Group of likely synonymous drugs',
};

sub consume {
    my $self = shift;
    my @groups = @_;
    for my $group (@groups){
        for my $bridge($group->bridges){
             $self->add_bridge(drug_id => $bridge->drug_id);
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
