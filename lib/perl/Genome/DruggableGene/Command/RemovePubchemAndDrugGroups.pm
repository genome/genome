package Genome::DruggableGene::Command::RemovePubchemAndDrugGroups;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::RemovePubchemAndDrugGroups {
    is => 'Genome::Command::Base',
    has => [
        version => {
            is => 'Text',
            is_input => 1,
            doc => 'Version identifier for the infile (ex 3)',
        },
    ],
    doc => 'Remove existing drug/pubchem groups'
};

sub help_brief {
      'Remove existing drug/pubchem groups'
}

sub execute {
    my $self = shift;
    Genome::DruggableGene::Command::DrugNameGroup::RemoveAll->execute(); #remove drug groups
    $self->remove_pubchem; #remove pubchem
    return 1;
}

sub remove_pubchem {
    my $self = shift;
    my @drugs = Genome::DruggableGene::DrugNameReport->get(source_db_name => 'PubChem', source_db_version => $self->version);
    for my $drug (@drugs){
        for my $alt_name ($drug->drug_alt_names){
            $alt_name->delete;
        }
        for my $cat ($drug->drug_categories){
            $cat->delete;
        }
        $drug->delete;
    }
    return 1;
}

1;
