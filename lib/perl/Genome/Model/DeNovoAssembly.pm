package Genome::Model::DeNovoAssembly;

use strict;
use warnings;

use Genome;
use Genome::ProcessingProfile::DeNovoAssembly;

class Genome::Model::DeNovoAssembly {
    is => 'Genome::ModelDeprecated',
    has => [
        center_name => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'center_name', value_class_name => 'UR::Value' ],
            is_mutable => 1,
            default_value => 'wugc',
            doc => 'Center name',
        },
        import_location => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'import_location', value_class_name => 'UR::Value' ],
            is_mutable => 1,
            is_optional => 1,
            doc => 'Directory to import assembly files from',
        },
	    map({
                $_ => {
                    via => 'processing_profile',
                }
            } Genome::ProcessingProfile::DeNovoAssembly->params_for_class
        ),
    ],
};

sub build_subclass_name {
    return 'de-novo-assembly';
}

sub _additional_parts_for_default_name {
    my $self = shift;

    my @parts;
    my $subject = $self->subject;
    if ( $subject->isa('Genome::Sample') and defined $subject->tissue_desc ) {
        my $tissue_name_part = $self->_get_name_part_from_tissue_desc($subject->tissue_desc);
        push @parts, $tissue_name_part if defined $tissue_name_part;
    }

    my $center_name = $self->center_name;
    if ( not $center_name ) {
        Carp::confess('No center name to get default model name for de novo assembly model.');
    }
    if ( $center_name ne 'WUGC' ) {
        push @parts, $center_name;
    }
    if ( $self->import_location ) {
        push @parts, $self->import_location;
    }

    return @parts;
}

sub _get_name_part_from_tissue_desc {
    my ($self, $tissue_desc) = @_;

    if ( $tissue_desc =~ /dna_/i ) {
        # zs5_g_dna_r_retroauricular crease
        my @tokens = split(/\_/, $tissue_desc);
        return if not @tokens;
        return Genome::Utility::Text::capitalize_words($tokens[$#tokens]);
    }

    return if $tissue_desc =~ /^\d/; # starts w/ number
    return if $tissue_desc =~ /^\w\d+$/; # starts w/ one letter then numbers

    $tissue_desc =~ s/,//g;
    return Genome::Utility::Text::capitalize_words($tissue_desc);
}

sub sequencing_platform { return; }

sub instrument_data_class_name { return 'Genome::InstrumentData'; }

sub compatible_instrument_data_subclasses {
    return (
        'Genome::InstrumentData::Solexa',
        'Genome::InstrumentData::Imported'
    );
}

sub compatible_instrument_data {
    my $self = shift;
    # assign only solexa and imported inst data
    my @all_data = $self->SUPER::compatible_instrument_data(@_);
    my @compatible_data;
    for my $data ( @all_data ) {
        next if not grep {$data->subclass_name eq $_} $self->compatible_instrument_data_subclasses;
        push @compatible_data, $data;
    }
    return @compatible_data;
}

1;

