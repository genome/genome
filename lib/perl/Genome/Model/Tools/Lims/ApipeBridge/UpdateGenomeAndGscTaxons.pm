package Genome::Model::Tools::Lims::ApipeBridge::UpdateGenomeAndGscTaxons;

use strict;
use warnings;

use Genome;

use Data::Dumper;

class Genome::Model::Tools::Lims::ApipeBridge::UpdateGenomeAndGscTaxons { 
    is => 'Command::V2',
    has => [
        taxons => {
            is => 'Genome::Taxon',
            is_many => 1,
            doc => 'Taxons to fix.',
            shell_args_position => 1,
        },
        domain => {
            is => 'Text',
            is_optional => 1,
            valid_values => _valid_values_for('domain'),
            doc => _doc_for('domain'),
        },
        estimated_genome_size => {
            is => 'Number',
            is_optional => 1,
            doc => _doc_for('estimated_genome_size'),
        },
        gram_stain_category => {
            is => 'Text',
            is_optional => 1,
            valid_values => _valid_values_for('gram_stain_category'),
            doc => _doc_for('gram_stain_category'),
        },
    ],
};

sub help_brief { return 'Update a property of a taxon in Genome and LIMS'; }
sub help_detail { return help_brief(); }

sub _property_for {
    my $property_name = shift;
    my $property = Genome::Taxon->__meta__->property_meta_for_name($property_name);
    Carp::confess('Failed to get property for name! '.$property_name) if not $property;
    return $property;
}

sub _doc_for {
    my $property_name = shift;
    my $property = _property_for($property_name);
    return $property->doc;
}

sub _valid_values_for {
    my $property_name = shift;
    my $property = _property_for($property_name);
    return $property->valid_values;
}

sub execute {
    my $self = shift;
    $self->debug_message('Update genome and gsc taxons...');

    my @property_names = (qw/ domain estimated_genome_size gram_stain_category /);
    my @property_names_to_update = grep { defined $self->$_ } @property_names;
    if ( not @property_names_to_update ) {
        $self->error_message('No properties to update!');
        return;
    }

    for my $taxon ( $self->taxons ) {
        $self->debug_message('Taxon: '.join(' ', map { $taxon->$_ } (qw/ id name /)));

        my $gsc_taxon = Genome::Site::TGI::Taxon->get(id => $taxon->id);
        if ( not $gsc_taxon ) {
            $self->error_message('Cannot update taxon! Failed to get GSC taxon for id: '.$taxon->id.'. Skipping...');
            next;
        }

        for my $name ( @property_names_to_update ) {
            my $value = $self->$name;
            $self->debug_message(
                sprintf(
                    " Update genome %s from %s to %s",
                    $name, ($taxon->$name || 'NULL'), $value,
                )
            );
            $taxon->$name($value) if not defined $taxon->$name or $taxon->$name ne $value;

            $self->debug_message(
                sprintf(" Update gsc %s from %s to %s\n",
                    $name, ($gsc_taxon->$name || 'NULL'), $value,
                )
            );
            $gsc_taxon->$name($value) if not defined $gsc_taxon->$name or $gsc_taxon->$name ne $value;
        }
    }

    $self->debug_message('Done');
    return 1;
}

1;

