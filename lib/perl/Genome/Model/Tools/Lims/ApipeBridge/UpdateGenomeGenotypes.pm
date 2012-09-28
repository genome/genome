package Genome::Model::Tools::Lims::ApipeBridge::UpdateGenomeGenotypes;

use strict;
use warnings;

class Genome::Model::Tools::Lims::ApipeBridge::UpdateGenomeGenotypes {
    is => 'Command::V2',
    has => [
        genotype_ids => {
            is => 'Text',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The id of the genotype to update.',
        },
    ],
};

sub help_brief { return 'Update APIPE genotype instrument data from LIMS'; }
sub help_detail { return help_brief(); }

sub execute {
    my $self = shift;
    $self->status_message('Updage genome genotype instrument data...');

    my ($attempted, $success) = (qw/ 0 0 /);
    for my $genotype_id ( $self->genotype_ids ) {
        $attempted++;
        $success++ if $self->_update_genotype($genotype_id);
    }

    $self->status_message('Attempted: '.$attempted);
    $self->status_message('Success: '.$success);
    $self->status_message('Done');
    return 1;
}

sub _update_genotype {
    my ($self, $genotype_id) = @_;

    my $instrument_data = Genome::InstrumentData->get($genotype_id);
    if ( not $instrument_data ) {
        $self->status_message('No inst data for id! '.$genotype_id);
        return;
    }

    my $genotype = Genome::Site::TGI::Synchronize::Classes->get($genotype_id);
    if ( not $genotype) {
        $self->error_message("No genotype for id! $genotype_id");
        return;
    }

    my $chip_name = $genotype->chip_name;
    if ( not $chip_name ) { # some platforms do not have this information
        $self->status_message( join(' ', $instrument_data->id, 'chip-type-na', 'na', 'version-na', 'na' ) );
        return 1;
    }
    my $version = $genotype->version;
    $version ||= 'na';

    my $chip_attribute = $instrument_data->attributes(attribute_label => 'chip_name');
    my $chip_status = 'ok';
    if ( not $chip_attribute ) {
        $chip_status = 'create';
        $chip_attribute = $instrument_data->add_attribute(attribute_label => 'chip_name', attribute_value => $chip_name);
    }
    elsif ( $chip_attribute->attribute_value ne $chip_name ) {
        $chip_status = 'update';
        $chip_attribute->delete;
        $chip_attribute = $instrument_data->add_attribute(attribute_label => 'chip_name', attribute_value => $chip_name);
    }
    if ( not $chip_attribute or $chip_attribute->attribute_value ne $chip_name ) {
        $self->error_message('Failed to add attribute for chip_name! '.$genotype_id);
        return;
    }

    my $version_attribute = $instrument_data->attributes(attribute_label => 'version');
    my $verion_status = 'ok';
    if ( not $version_attribute ) {
        $verion_status = 'create';
        $version_attribute = $instrument_data->add_attribute(attribute_label => 'version', attribute_value => $version);
    }
    elsif ( $version_attribute->attribute_value ne $version ) {
        $verion_status = 'update';
        $version_attribute->delete;
        $version_attribute = $instrument_data->add_attribute(attribute_label => 'version', attribute_value => $version);
    }
    if ( not $version_attribute or $version_attribute->attribute_value ne $version ) {
        $self->error_message('Failed to add attribute for version! '.$genotype_id);
        return;
    }

    $self->status_message( 
        join(' ', $instrument_data->id, 'chip-'.$chip_status, $chip_attribute->attribute_value, 'version-'.$verion_status, $version_attribute->attribute_value)
    );

    return 1;
}

1;

