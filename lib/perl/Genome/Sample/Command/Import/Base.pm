package Genome::Sample::Command::Import::Base;

use strict;
use warnings;

use Genome;

class Genome::Sample::Command::Import::Base {
    is => 'Command',
    is_abstract => 1,
    doc => 'Import samples from known sources',
    has => [
        _taxon => { is => 'Genome::Taxon', is_optional => 1, },
        _individual => { is => 'Genome::Individual', is_optional => 1, },
        _sample => { is => 'Genome::Sample', is_optional => 1, },
        _library => { is => 'Genome::Library', is_optional => 1, },
        _created_objects => { is => 'ARRAY', is_optional => 1, },
    ],
};

sub help_brief {
    return 'import samples from known sources';
}

sub help_detail {
    return help_brief();
}

sub _update_object {
    my ($self, $obj, %params) = @_;

    $self->status_message('Update '.$obj->name.' ('.$obj->id.')');
    my $force = delete $params{__force__};
    $self->status_message('Force is '.($force ? 'on' : 'off'));
    $self->status_message('Params: '._display_string_for_params(\%params));

    for my $attr ( keys %params ) {
        my $val = $obj->$attr;
        if ( defined $val and not $force ) {
            $self->status_message("Not updating '$attr' for ".$obj->id." because it already has a value ($val)");
            next;
        }
        $obj->$attr( $params{$attr} );
    }

    $self->status_message('Update...OK');
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
    my $library_ext = delete $params{library};
    Carp::confess('No library extention given to import') if not $library_ext;

    # taxon
    $self->_taxon( Genome::Taxon->get(name => $taxon_name) );
    Carp::confess("Cannot get taxon for '$taxon_name'") if not $self->_taxon;
    $self->status_message('Found taxon: '.$self->_taxon->__display_name__);

    # sample
    my $sample = Genome::Sample->get(name => $sample_name);
    if ( $sample ) {
        $self->_sample($sample);
        $self->status_message('Found sample: '.join(' ', map{ $sample->$_ } (qw/ id name/)));
        if ( %$sample_params ) { # got additional params - try to update
            my $update = $self->_update_object($sample, %$sample_params);
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
        $self->_bail('Sample ('.$sample->id.') source id ('.$sample->source_id.') does not match found individual ('.$individual->id.')');
        return;
    }

    if ( not $sample->source_type ) {
        $sample->source_type( $individual->subject_type );
    }
    if ( $sample->source_id ne $individual->id ) {
        $self->_bail('Sample ('.$sample->id.') source type ('.$sample->source_type.') does not match individual ('.$individual->subject_type.')');
        return;
    }

    # library
    my $library = $self->_get_or_create_library_for_extension($library_ext);
    return if not $library;

    return 1;
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
    my $min_unique_name_parts = $self->_minimum_unique_name_parts - 1;
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
            $self->_bail("Found multiple individuals for similar samples: ".join(' ', keys %individuals));
            return;
        }
        my ($individual) = values %individuals;
        return $self->_individual($individual);
    }

    my $individual_for_given_upn = Genome::Individual->get(upn => $upn);
    return if not $individual_for_given_upn;
    return $self->_individual($individual_for_given_upn);
}

sub _minimum_unique_name_parts {
    return 2;
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
        $self->_bail('Could not create individual');
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
        $self->_bail('Cannot create sample');
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
        $self->_bail('Cannot not create library to import sample');
        return;
    }

    my $created_objects = $self->_created_objects;
    push @$created_objects, $library;
    $self->_created_objects($created_objects);

    $self->status_message('Create library: '.join(' ', map{ $library->$_ } (qw/ id name/)));
    
    return $self->_library($library);
}

sub _bail {
    my ($self, $msg) = @_;

    $self->error_message($msg);

    my $created_objects = $self->_created_objects;
    return if not defined $created_objects;

    $self->status_message('Encountered an error, removing newly created objects.');

    for my $obj ( @$created_objects ) { 
        my $transaction = UR::Context::Transaction->begin();
        $obj->delete;
    }

    return 1;
}

sub _display_string_for_params {
    my ($params) = shift;

    my $string = "\n";
    for my $key ( sort keys %$params ) {
        $string .= sprintf(" %s: %s\n", $key, ( ref($params->{$key}) ? $params->{$key}->__display_name__ : $params->{$key} ));
    }

    return $string;
}

1;

