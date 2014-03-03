package Genome::Sample::Command::Import::Base;

use strict;
use warnings;

use Genome;

class Genome::Sample::Command::Import::Base {
    is => 'Command',
    is_abstract => 1,
    doc => 'Import samples from known sources',
    has => [
        name => {
            is => 'Text',
            doc => 'Sample name. Source name will be derived from the sample name. The typical format is "NOMENCLATURE-SOURCE_ID/NAME-SAMPLE_ID/NAME, but may vary based on source.',
        },
        extraction_type => {
            is => 'Text',
            default_value => 'genomic dna',
            valid_values => [ 'cdna', 'genomic dna', 'ipr product', 'rna', 'total rna', ],
            doc => 'The extraction type of the sample.',
        },
    ],
    has_optional => [
        individual_attributes => {
            is => 'Text',
            is_many => 1,
            doc => 'Additional attributes to add to the individual. Give as key value pairs. Separate key and value with an equals (=) and pairs with a comma (,). Ex: attr1=val1,attr2=val2',
        },
        sample_attributes => {
            is => 'Text',
            is_many => 1,
            doc => 'Additional attributes to add to the sample. Give as key value pairs. Separate key and value with an equals (=) and pairs with a comma (,). Ex: attr1=val1,attr2=val2',
        },
        library_attributes => {
            is => 'Text',
            is_many => 1,
            doc => 'Additional attributes to add to the library. Give as key value pairs. Separate key and value with an equals (=) and pairs with a comma (,). Ex: attr1=val1,attr2=val2',
        },
    ],
    has_optional_transient => [
        # taxon
        _taxon => { is => 'Genome::Taxon', is_optional => 1, },
        # source
        _individual => { is => 'Genome::Individual', is_optional => 1, },
        _individual_name => { is => 'Text', },
        # sample
        _sample => { is => 'Genome::Sample', is_optional => 1, },
        # library
        _library => { is => 'Genome::Library', is_optional => 1, },
        # misc
        _created_objects => { is => 'ARRAY', is_optional => 1, },
        _minimum_unique_source_name_parts => { is => 'Number', default_value => 2, },
    ],
};

sub help_detail {
    return;
}

sub execute {
    my $self = shift;
    $self->status_message('Import '.$self->nomenclature.' sample...');

    my $individual_name_ok = $self->_validate_name_and_set_individual_name;
    return if not $individual_name_ok;

    my %individual_attributes = $self->_resolve_individual_attributes;
    return if not %individual_attributes;

    my %sample_attributes = $self->_resolve_sample_attributes;
    return if not %sample_attributes;
    
    my %library_attributes = $self->_resolve_library_attributes;

    my $import = $self->_import(
        taxon => $self->taxon_name,
        individual => \%individual_attributes,
        sample => \%sample_attributes,
        library => \%library_attributes,
    );
    return if not $import;

    $self->status_message('Import sample...OK');
    return 1;
}

sub _import {
    my ($self, %params) = @_;

    # params
    Carp::confess('No params given to import') if not %params;
    my $taxon_name = delete $params{taxon};
    Carp::confess('No taxon name given to import') if not $taxon_name;
    my $individual_params = delete $params{individual};
    Carp::confess('No individual params given to import') if not $individual_params;
    my $individual_upn = delete $individual_params->{upn};
    Carp::confess('No individual upn in individual params given to import') if not $individual_upn;
    my $sample_params = delete $params{sample};
    Carp::confess('No sample params given to import') if not $sample_params;
    my $sample_name = delete $sample_params->{name};
    Carp::confess('No sample name in sample params given to import') if not $sample_name;
    my $library_params = delete $params{library};
    Carp::confess('No library params given to import') if not $library_params;
    my $library_ext = delete $library_params->{ext};
    Carp::confess('No library ext given to import') if not $library_ext;

    # taxon
    $self->_taxon( Genome::Taxon->get(name => $taxon_name) );
    Carp::confess("Cannot get taxon for '$taxon_name'") if not $self->_taxon;
    $self->status_message('Found taxon: '.$self->_taxon->__display_name__);

    # sample
    my $sample = Genome::Sample->get(name => $sample_name);
    if ( $sample ) {
        $self->_sample($sample);
        $self->status_message('Found sample: '.join(' ', map{ $sample->$_ } (qw/ id name/)));
        if ( %$sample_params ) { # got additional attributes - try to update
            my $update = $self->_update_attributes($sample, %$sample_params);
            return if not $update;
        }
    }
    else { # create, set individual later
        $sample = $self->_create_sample(
            name => $sample_name,
            %$sample_params,
        );
        return if not $sample;
    }

    # individual
    my $individual = $self->_get_individual($individual_upn); # get by sample and upn
    if ( not $individual ) {
        $individual = $self->_create_individual(
            upn => $individual_upn,
            %$individual_params,
        );
        return if not $individual;
    }
    else {
        $self->status_message('Found individual: '.join(' ', map{ $individual->$_ } (qw/ id name upn /)));
        $self->status_message("Found individual 'upn' (".$self->_individual->upn.") does not match the upn given ($individual_upn). This is probably ok.") if $individual->upn ne $individual_upn;
    }

    if ( not $sample->source_id ) {
        $sample->source_id( $individual->id );
    }
    if ( $sample->source_id ne $individual->id ) {
        $self->error_message('Sample ('.$sample->id.') source id ('.$sample->source_id.') does not match found individual ('.$individual->id.')');
        return;
    }

    if ( not $sample->source_type ) {
        $sample->source_type( $individual->subject_type );
    }
    if ( $sample->source_id ne $individual->id ) {
        $self->error_message('Sample ('.$sample->id.') source type ('.$sample->source_type.') does not match individual ('.$individual->subject_type.')');
        return;
    }

    # library
    my $library = $self->_get_or_create_library_for_extension($library_ext);
    $library = $self->_set_library_params($library, $library_params);
    return if not $library;

    return 1;
}

sub _set_library_params {
    my $self = shift;
    my $library = shift;
    my $params = shift;

    for my $param_name (keys %{$params}) {
        $library->$param_name($params->{$param_name});
    }
    return $library;
}

sub _validate_name_and_set_individual_name {
    my $self = shift;
    $self->status_message('Validate sample name and resolve individual name...');

    my $name_regexp = $self->name_regexp;
    my $name = $self->name;
    $self->status_message('Sample name: '.$name);
    $self->status_message('Sample regexp: '.$name_regexp);
    if ( $name !~ /$name_regexp/ ) {
        $self->error_message("Name ($name) is invalid!");
        return;
    }
    my $individual_name = $1;
    $self->status_message('Individual name: '.$individual_name);
    $self->_individual_name($individual_name);

    $self->status_message('Validate sample name and resolve individual name...done');
    return 1
}

sub _get_individual {
    my ($self, $upn) = @_;

    my $sample = $self->_sample;
    Carp::confess('No sample set to get individual') if not $sample;
    Carp::confess('No "upn" given to get individual') if not $upn;

    if ( my $individual = $sample->source ) {
        return $self->_individual($individual);
    }

    my %individuals_from_similar_samples;
    my @tokens = split('-', $sample->name);
    my $min_unique_name_parts = $self->_minimum_unique_source_name_parts - 1;
    for ( my $i = $#tokens - 1; $i > $min_unique_name_parts; $i--  ) { # go down to 2 levels
        my $extraction_label = join('-', @tokens[0..$i]);
        my @samples = Genome::Sample->get('extraction_label like' => $extraction_label.'%');
        for my $sample ( @samples ) {
            my $source = $sample->source;
            next if not $source;
            next if not $source->isa('Genome::Individual');
            $individuals_from_similar_samples{$source->id} = $source;
        }
        last if %individuals_from_similar_samples;
    }

    if ( %individuals_from_similar_samples ) {
        my %individuals = map { $_->id => $_ } values %individuals_from_similar_samples;
        if ( keys %individuals > 1 ) {
            $self->error_message("Found multiple individuals for similar samples: ".join(' ', keys %individuals));
            return;
        }
        my ($individual) = values %individuals;
        return $self->_individual($individual);
    }

    my $individual_for_given_upn = Genome::Individual->get(upn => $upn);
    return if not $individual_for_given_upn;
    return $self->_individual($individual_for_given_upn);
}

sub _create_individual {
    my ($self, %params) = @_;

    Carp::confess('No "upn" given to create individual') if not $params{upn};
    Carp::confess('No "nomenclature" given to create individual') if not $params{nomenclature};
    Carp::confess('No taxon set to create individual') if not $self->_taxon;

    $params{name} = $params{upn} if not $params{name};
    $params{taxon} = $self->_taxon;
    $params{gender} = 'unspecified' if not $params{gender};

    $self->status_message('Create individual...');
    $self->status_message('Individual params: '._display_string_for_params(\%params));
    my $transaction = UR::Context::Transaction->begin();
    my $individual = Genome::Individual->create(%params);
    if ( not defined $individual ) {
        $self->error_message('Could not create individual');
        return;
    }

    my $created_objects = $self->_created_objects;
    push @$created_objects, $individual;
    $self->_created_objects($created_objects);
    $self->status_message('Create individual: '.join(' ', map{ $individual->$_ } (qw/ id name/)));

    return $self->_individual($individual);
}

sub _create_sample {
    my ($self, %params) = @_;

    Carp::confess('No name given to create sample') if not $params{name};
    Carp::confess('No nomenclature set to create sample') if not $params{nomenclature};

    if ( $self->_individual ) {
        $params{source_id} = $self->_individual->id;
        $params{source_type} = $self->_individual->subject_type;
    }

    $self->status_message('Create sample...');
    $self->status_message('Sample params: '._display_string_for_params(\%params));
    my $sample = Genome::Sample->create(%params);
    if ( not defined $sample ) {
        $self->error_message('Cannot create sample');
        return;
    }

    my $created_objects = $self->_created_objects;
    push @$created_objects, $sample;
    $self->_created_objects($created_objects);

    $self->status_message('Create sample: '.join(' ', map { $sample->$_ } (qw/ id name /)));

    return $self->_sample($sample);
}

sub _get_or_create_library_for_extension {
    my ($self, $ext) = @_;

    my $library = $self->_get_library_for_extension($ext);
    return $library if $library;

    return $self->_create_library_for_extension($ext);
}

sub _get_library_name_for_extension {
    my ($self, $ext) = @_;

    Carp::confess('No sample set to get or create library') if not $self->_sample;
    Carp::confess('No library extension') if not defined $ext;
    my @valid_exts = (qw/ extlibs microarraylib /);
    Carp::confess("Invalid library extension ($ext). Valid extentions: ".join(' ', @valid_exts)) if not grep { $ext eq $_ } @valid_exts;

    return $self->_sample->name.'-'.$ext;
}

sub _get_library_for_extension {
    my ($self, $ext) = @_;

    my $name = $self->_get_library_name_for_extension($ext); # confess on error
    my $library = Genome::Library->get(name => $name);
    return if not $library;

    $self->status_message('Found library: '.join(' ', map{ $library->$_ } (qw/ id name/)));

    return $self->_library($library);

}

sub _create_library_for_extension {
    my ($self, $ext) = @_;

    my %params = (
        name => $self->_get_library_name_for_extension($ext), # confess on error
        sample_id => $self->_sample->id,
    );

    $self->status_message('Create library...');
    $self->status_message('Library params: '._display_string_for_params(\%params));
    my $library = Genome::Library->create(%params);
    if ( not $library ) {
        $self->error_message('Cannot not create library to import sample');
        return;
    }

    my $created_objects = $self->_created_objects;
    push @$created_objects, $library;
    $self->_created_objects($created_objects);

    $self->status_message('Create library: '.join(' ', map{ $library->$_ } (qw/ id name/)));
    
    return $self->_library($library);
}

sub _display_string_for_params {
    my ($params) = shift;

    my $string = "\n";
    for my $key ( sort keys %$params ) {
        $string .= sprintf(" %s: %s\n", $key, ( ref($params->{$key}) ? $params->{$key}->__display_name__ : $params->{$key} ));
    }

    return $string;
}

sub _resolve_individual_attributes {
    my $self = shift;
    my %attributes = (
        nomenclature => $self->nomenclature,
        name => $self->_individual_name,
        upn => $self->_individual_name,
    );
    return if not $self->_resolve_attributes('individual', \%attributes);
    return %attributes;
}

sub _resolve_sample_attributes {
    my $self = shift;
    my %attributes = (
        nomenclature => $self->nomenclature,
        name => $self->name,
        extraction_label => $self->name,
        cell_type => 'unknown',
    );
    return if not $self->_resolve_attributes('sample', \%attributes);
    $attributes{extraction_type} = $self->extraction_type if not defined $attributes{extraction_type};
    return %attributes;
}

sub _resolve_library_attributes {
    my $self = shift;
    my %attributes = (
        ext => "extlibs",
    );
    return if not $self->_resolve_attributes('library', \%attributes);
    return %attributes;
}

sub _resolve_attributes {
    my ($self, $type, $attributes) = @_;

    my $attributes_method = $type.'_attributes';
    my @additional_attributes = $self->$attributes_method;
    if ( @additional_attributes ) {
        no warnings;
        my %additional_attributes = map { split('=') } @additional_attributes;
        use warnings;

        for my $label ( keys %additional_attributes ) {
            if ( defined $additional_attributes{$label} ) {
                $attributes->{$label} = $additional_attributes{$label};
            }
            else {
                $self->error_message(ucfirst($type)." attribute label ($label) does not have a value!");
                return;
            }
        }
    }

    my $attribute_names_method = '_'.$type.'_attribute_names';
    my $names = eval{ $self->$attribute_names_method; };
    for my $name ( @$names ) {
        my $value = eval{ $self->$name; };
        if ( $@ ) {
            $self->error_message($@);
            return;
        }
        next if not defined $value;
        $name =~ s/^$type\_//; # attributes may have the same name, like common_name, so remove the type in front
        $attributes->{$name} = $value;
    }

    return 1;
}

sub _update_attributes {
    my ($self, $obj, %attributes) = @_;

    $self->status_message('Update '.$obj->name.' ('.$obj->id.')');
    my $force = delete $attributes{__force__};
    $self->status_message('Force is '.($force ? 'on' : 'off'));
    $self->status_message('Params: '._display_string_for_params(\%attributes));

    for my $label ( keys %attributes ) {
        my $value = eval{ $obj->attributes(attribute_label => $label)->attribute_value; };
        if ( defined $value ) {
            $self->status_message("Not updating '$label' for ".$obj->id." because it already has a value ($value)");
            next;
        }
        $obj->add_attribute(
            attribute_label => $label,
            attribute_value => $attributes{$label},
            nomenclature => $self->nomenclature,
        );
    }

    $self->status_message('Update...OK');
    return 1;
}

1;

