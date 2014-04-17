package Genome::ModelDeprecated;

use strict;
use warnings;
use Genome;
use Regexp::Common;
use YAML;
use File::Path;

class Genome::ModelDeprecated {
    is => 'Genome::Model',
    is_abstract => 1,
    has => [
        subject_name => {
            via => 'subject',
            to => 'name',
        },
        _sample_subject => {
            is => 'Genome::Sample',
            is_optional => 1,
            id_by => 'subject_id',
        },
        subject_type => { 
            is => 'Text', 
            valid_values => ["species_name","sample_group","sample_name"], 
            calculate_from => 'subject_class_name',
            calculate => q|
                #This could potentially live someplace else like the previous giant hash
                my %types = (
                    'Genome::Sample' => 'sample_name',
                    'Genome::PopulationGroup' => 'sample_group',
                    'Genome::Individual' => 'sample_group',
                    'Genome::Taxon' => 'species_name',
                );
                return $types{$subject_class_name};
            |, 
        },
        processing_profile_name => { via => 'processing_profile', to => 'name' },
    ],
    has_optional => [
        is_default => { 
            is => 'Boolean',
            doc => 'flag the model as the default system "answer" for its subject'
        },
#        auto_assign_inst_data => {
#            is => 'Boolean',
#            calculate_from => ['_auto_assign_inst_data'],
#            calculate => q{ $_auto_assign_inst_data; }
#        },
#        auto_build_alignments => {
#            is => 'Boolean',
#            calculate_from => ['_auto_build_alignments'],
#            calculate => q{ $_auto_build_alignments; }
#        }, # TODO: rename to auto_build
#        apipe_cron_status => {
#            via => 'notes',
#            to => 'body_text',
#            where => [ header_text => 'apipe_cron_status' ],
#            is_mutable => 0,
#        },
    ],
    has_optional_many => [
        group_ids => { via => 'model_groups', to => 'id' },
        group_names => { via => 'model_groups', to => 'name' },
        project_names => { is => 'Text', via => 'projects', to => 'name', },

        # model links are now deprecated in favor of making the from_model an
        # input on the to_model
        from_model_links => { 
            is => 'Genome::Model::Link', reverse_as => 'to_model',
            doc => 'bridge table entries where this is the "to" model (used to retrieve models this model is "from"' 
        },
        from_models => { 
            is => 'Genome::Model', via => 'from_model_links', to => 'from_model',
            doc => 'Genome models that contribute "to" this model' 
        },
        to_model_links => { 
            is => 'Genome::Model::Link', reverse_as => 'from_model',
            doc => 'bridge entries where this is the "from" model(used to retrieve models models this model is "to")' 
        },
        to_models => { 
            is => 'Genome::Model', via => 'to_model_links', to => 'to_model',
            doc => 'Genome models this model contributes "to"' 
        },
        attributes => { 
            is => 'Genome::MiscAttribute', 
            reverse_as => '_model', 
            where => [ entity_class_name => 'Genome::Model' ] 
        },
        # TODO: these go into a model subclass for models which apply directly to sequencer data
        sequencing_platform         => { via => 'processing_profile' },
        instrument_data_class_name  => { is => 'Text',
                                        calculate_from => 'sequencing_platform',
                                        calculate => q| 'Genome::InstrumentData::' . ucfirst($sequencing_platform) |,
                                        doc => 'the class of instrument data assignable to this model' },
        instrument_data_inputs => {
            is => 'Genome::Model::Input',
            reverse_as => 'model',
            where => [ name => 'instrument_data' ],
        },
    ],
    has_optional_deprecated => [
        events => { 
            # TODO: events are used in old-style processing profiles for workflow steps, 
            # which is deprecated, but also have a non-deprecated uses like logging
            # model creation, addition of instrument data or other input changes, etc.,
            # abandon dates, and other things we don't want to give special fields to.
            is => 'Genome::Model::Event', 
            reverse_as => 'model', 
            is_many => 1,
            doc => 'all events which have occurred for this model' 
        },

        # this is all old junk but things really use them right now 
        reports                 => { via => 'last_succeeded_build' },
        reports_directory       => { via => 'last_succeeded_build' },

        # these go on refalign models
        region_of_interest_set_name => { 
            is => 'Text',
            is_many => 1, 
            is_mutable => 1,
            via => 'inputs', 
            to => 'value_id',
            where => [ name => 'region_of_interest_set_name', value_class_name => 'UR::Value' ], 
        },
        merge_roi_set => {
            is_mutable => 1,
            via => 'inputs', 
            to => 'value_id',
            where => [ name => 'merge_roi_set', value_class_name => 'UR::Value' ], 
        },
        short_roi_names => {
            is_mutable => 1,
            via => 'inputs', 
            to => 'value_id',
            where => [ name => 'short_roi_names', value_class_name => 'UR::Value' ], 
        },
        roi_track_name => {
            is_mutable => 1,
            via => 'inputs', 
            to => 'value_id',
            where => [ name => 'roi_track_name', value_class_name => 'UR::Value' ], 
        },
    ],
    has_optional_calculated => [
        individual_common_name => {
            is => 'Text',
            calculate_from => 'subject',
            calculate => q{
                if ($subject->class eq 'Genome::Individual') {
                    return $subject->common_name();
                } elsif($subject->class eq 'Genome::Sample') {
                    return $subject->patient_common_name();
                } else {
                    return undef;
                }
            },
        },
        sample_names => {
            is => 'Array',
            calculate => q{
                my @s = $self->get_all_possible_samples();
                return sort map {$_->name()} @s;
            },
        },
    ],
    has_deprecated_optional => [
        last_complete_build_directory    => { calculate => q|$b = $self->last_complete_build; return unless $b; return $b->data_directory| },
        last_succeeded_build_directory   => { calculate => q|$b = $self->last_succeeded_build; return unless $b; return $b->data_directory| },
        build_statuses                  => { via => 'builds', to => 'master_event_status', is_many => 1 },
        build_ids                       => { via => 'builds', to => 'id', is_many => 1 },
        sample_name                     => { is => 'Text', doc => 'deprecated column with explicit sample_name tracking' },
    ],
};

sub create {
    my $class = shift;
    my $bx = $class->define_boolexpr(@_);
    if ($bx->specifies_value_for("subject_class_name")) {
        $bx = $bx->remove_filter("subject_class_name");
    }
    my $self = $class->SUPER::create($bx);
    #if (my $subject = $self->subject) {
    #    $self->subject_class_name($subject->class);
    #}
    return $self;
}

# TODO This can likely be simplified once all subjects are a subclass of Genome::Subject
#If a user defines a model with a name (and possibly type), we need to find/make sure there's an
#appropriate subject to use based upon that name/type.
sub _resolve_subject_from_name_and_type {
    my $class = shift;
    my $subject_name = shift;
    my $subject_type = shift;

    if (not defined $subject_name) {
        $class->error_message("bad data--missing subject_name!");
        return;
    }

    my $try_all_types = 0;
    if (not defined $subject_type) {
        #If they didn't give a subject type, we'll keep trying subjects until we find something that sticks.
        $try_all_types = 1;
    }

    my @subjects = ();

    if($try_all_types or $subject_type eq 'sample_name') {
        my $subject = Genome::Sample->get(name => $subject_name);
        return $subject if $subject; #sample_name is the favoured default.  If we get one, use it.
    }
    if ($try_all_types or $subject_type eq 'species_name') {
        push @subjects, Genome::Taxon->get(name => $subject_name);
    }
    if ($try_all_types or $subject_type eq 'library_name') {
        push @subjects, Genome::Library->get(name => $subject_name);
    }
    if ($try_all_types or $subject_type eq 'genomic_dna') {
        push @subjects, Genome::Sample->get(extraction_label => $subject_name, extraction_type => 'genomic dna');
    }

    #This case will only be entered if the user asked specifically for a sample_group
    if ($subject_type and $subject_type eq 'sample_group') {
        push @subjects,
            Genome::Individual->get(name => $subject_name),
            Genome::ModelGroup->get(name => $subject_name),
            Genome::PopulationGroup->get(name => $subject_name);
    }

    if(scalar @subjects == 1) {
        return $subjects[0];
    } elsif (scalar @subjects) {
        my $null = '<NULL>';
        $class->error_message('Multiple matches for ' . join(', ',
            'subject_name: ' . ($subject_name || $null),
            'subject_type: ' . ($subject_type || $null),
        ) . '. Please specify a subject_type or use subject_id/subject_class_name instead.'
        );
        $class->error_message('Possible subjects named "' . $subject_name . '": ' . join(', ',
            map($_->class . ' #' . $_->id, @subjects)
        ));
    } else {
        #If we get here, nothing matched.
        my $null = '<NULL>';
        $class->error_message('Unable to determine a subject given ' . join(', ',
            'subject_name: ' . ($subject_name || $null),
            'subject_type: ' . ($subject_type || $null),
        ));
        return;
    }
}



sub get_all_possible_samples {
    my $self = shift;

    my @samples;
    if ( $self->subject_class_name eq 'Genome::Taxon' ) {
        my $taxon = Genome::Taxon->get(name => $self->subject_name);
        @samples = $taxon->samples();

        #data tracking is incomplete, so sometimes these need to be looked up via the sources
        my @sources = ($taxon->individuals, $taxon->population_groups);
        push @samples,
            map($_->samples, @sources);
    } elsif ($self->subject_class_name eq 'Genome::Sample'){
        @samples = ( $self->subject );
    } elsif ($self->subject_class_name eq 'Genome::Individual') {
        @samples = $self->subject->samples();
    #} elsif () {
        #TODO Possibly fill in for possibly Genome::PopulationGroup and possibly others (possibly)
    } else {
        @samples = ();
    }

    return @samples;
}

#< Instrument Data >#
sub input_for_instrument_data_id {
    my ($self, $id) = @_;
    return unless $id;
    for my $input ($self->instrument_data_inputs) {
        return $input if $input->value_id eq $id;
    }
    return;
}

sub input_for_instrument_data {
    my ($self, $instrument_data) = @_;
    return unless $instrument_data;
    return $self->input_for_instrument_data_id($instrument_data->id);
}

sub last_succeeded_build_id { return $_[0]->last_complete_build_id; }
sub last_complete_build_id {
    my $self = shift;
    my $last_complete_build = $self->last_complete_build;
    return unless $last_complete_build;
    return $last_complete_build->id;
}

sub abandoned_builds {
    my $self = shift;
    my @abandoned_builds = $self->builds_with_status('Abandoned');
    return @abandoned_builds;
}
sub failed_builds {
    my $self = shift;
    my @failed_builds = $self->builds_with_status('Failed');
    return @failed_builds;
}
sub running_builds {
    my $self = shift;
    my @running_builds = $self->builds_with_status('Running');
    return @running_builds;
}
sub scheduled_builds {
    my $self = shift;
    my @scheduled_builds = $self->builds_with_status('Scheduled');
    return @scheduled_builds;
}

sub current_running_build {
    my $self = shift;
    my @running_builds = $self->running_builds;
    my $current_running_build = pop(@running_builds);
    return $current_running_build;
}

sub current_running_build_id {
    my $self = shift;
    my $current_running_build = $self->current_running_build;
    unless ($current_running_build) {
        return;
    }
    return $current_running_build->id;
}

sub latest_build_directory {
    my $self = shift;
    my $current_running_build = $self->current_running_build;
    if (defined $current_running_build) {
        return $current_running_build->data_directory;
    }
    my $last_complete_build = $self->last_complete_build;
    if (defined $last_complete_build) {
        return $last_complete_build->data_directory;
    } else {
       return;
    }
}

sub yaml_string {
    my $self = shift;
    my $string = YAML::Dump($self);
    my @objects = $self->get_all_objects;
    for my $object (@objects) {
        $string .= $object->yaml_string;
    }
    return $string;
}

# TODO Will be removed when model links are phased out
# TODO please rename this -ss
sub add_to_model{
    my $self = shift;
    my (%params) = @_;
    my $model = delete $params{to_model};
    my $role = delete $params{role};
    $role||='member';

    $self->error_message("no to_model provided!") and die unless $model;
    my $from_id = $self->id;
    my $to_id = $model->id;
    unless( $to_id and $from_id){
        $self->error_message ( "no value for this model(from_model) id: <$from_id> or to_model id: <$to_id>");
        die;
    }
    my $reverse_bridge = Genome::Model::Link->get(from_model_id => $to_id, to_model_id => $from_id);
    if ($reverse_bridge){
        my $string =  "A model link already exists for these two models, and in the opposite direction than you specified:\n";
        $string .= "to_model: ".$reverse_bridge->to_model." (this model)\n";
        $string .= "from_model: ".$reverse_bridge->from_model." (the model you are trying to set as a 'to' model for this one)\n";
        $string .= "role: ".$reverse_bridge->role;
        $self->error_message($string);
        die;
    }
    my $bridge = Genome::Model::Link->get(from_model_id => $from_id, to_model_id => $to_id);
    if ($bridge){
        my $string =  "A model link already exists for these two models:\n";
        $string .= "to_model: ".$bridge->to_model." (the model you are trying to set as a 'to' model for this one)\n";
        $string .= "from_model: ".$bridge->from_model." (this model)\n";
        $string .= "role: ".$bridge->role;
        $self->error_message($string);
        die;
    }
    $bridge = Genome::Model::Link->create(from_model_id => $from_id, to_model_id => $to_id, role => $role);
    return $bridge;
}

# TODO Will be removed when model links are phased out
# TODO please rename this -ss
sub add_from_model{
    my $self = shift;
    my (%params) = @_;
    my $model = delete $params{from_model};
    my $role = delete $params{role};
    $role||='member';

    $self->error_message("no from_model provided!") and die unless $model;
    my $to_id = $self->id;
    my $from_id = $model->id;
    unless( $to_id and $from_id){
        $self->error_message ( "no value for this model(to_model) id: <$to_id> or from_model id: <$from_id>");
        die;
    }
    my $reverse_bridge = Genome::Model::Link->get(from_model_id => $to_id, to_model_id => $from_id);
    if ($reverse_bridge){
        my $string =  "A model link already exists for these two models, and in the opposite direction than you specified:\n";
        $string .= "to_model: ".$reverse_bridge->to_model." (the model you are trying to set as a 'from' model for this one)\n";
        $string .= "from_model: ".$reverse_bridge->from_model." (this model)\n";
        $string .= "role: ".$reverse_bridge->role;
        $self->error_message($string);
        die;
    }
    my $bridge = Genome::Model::Link->get(from_model_id => $from_id, to_model_id => $to_id);
    if ($bridge){
        my $string =  "A model link already exists for these two models:\n";
        $string .= "to_model: ".$bridge->to_model." (this model)\n";
        $string .= "from_model: ".$bridge->from_model." (the model you are trying to set as a 'from' model for this one)\n";
        $string .= "role: ".$bridge->role;
        $self->error_message($string);
        die;
    }
    $bridge = Genome::Model::Link->create(from_model_id => $from_id, to_model_id => $to_id, role => $role);
    return $bridge;
}
# This is called by both of the above.
sub _resolve_subclass_name_for_type_name {
    my ($class,$type_name) = @_;
    my @type_parts = split(' ',$type_name);

    my @sub_parts = map { ucfirst } @type_parts;
    my $subclass = join('',@sub_parts);

    my $class_name = join('::', 'Genome::Model' , $subclass);
    return $class_name;
}

sub _resolve_type_name_for_subclass_name {
    my ($class,$subclass_name) = @_;
    my ($ext) = ($subclass_name =~ /Genome::Model::(.*)/);
    return unless ($ext);
    my @words = $ext =~ /[a-z]+|[A-Z](?:[A-Z]+|[a-z]*)(?=$|[A-Z])/g;
    my $type_name = lc(join(" ", @words));
    return $type_name;
}

# TODO: please rename this -ss
sub get_all_objects {
    my $self = shift;
    my $sorter = sub { # not sure why we sort, but I put it in a anon sub for convenience
        return unless @_;
        if ( $_[0]->id =~ /^\-/) {
            return sort {$b->id cmp $a->id} @_;
        }
        else {
            return sort {$a->id cmp $b->id} @_;
        }
    };

    return map { $sorter->( $self->$_ ) } (qw{ inputs builds to_model_links from_model_links });
}

sub create_rule_limiting_instrument_data {
    my ($self, @instrument_data) = @_;
    @instrument_data = $self->instrument_data unless @instrument_data;
    return unless @instrument_data;

    # Find the smallest scale domain object that encompasses all the instrument data
    # and create a boolean expression for it.
    for my $accessor (qw/ library_id sample_id sample_source_id taxon_id /) {
        my @ids = map { $_->$accessor } @instrument_data;
        next unless @ids;
        next if grep { $_ ne $ids[0] } @ids;

        my $rule = $instrument_data[0]->define_boolexpr($accessor => $ids[0]);
        return $rule;
    }

    return;
}

sub latest_build_request_note {
    my $self = shift;
    my @notes = sort { $b->entry_date cmp $a->entry_date } grep { $_->header_text eq 'build_requested' } $self->notes;
    return unless @notes;
    return $notes[0];
}
    
sub time_of_last_build_request {
    my $self = shift;
    my $note = $self->latest_build_request_note;
    return unless $note;
    return $note->entry_date;
}

sub property_names_for_copy {
    my $class = shift;

    my $meta = eval{ $class->__meta__; };
    if ( not $meta ) {
        $class->error_message('Failed to get class meta for '.$class);
        return;
    }

    my @base_properties = (qw/
        auto_assign_inst_data auto_build_alignments processing_profile subject 
        /);

    my @input_properties = map { 
        $_->property_name
    } grep { 
        defined $_->via and $_->via eq 'inputs'
    } $meta->property_metas;

    return sort { $a cmp $b } ( @base_properties, @input_properties );
}

sub delete {
    my $self = shift;

    for my $group ($self->model_groups) {
        $self->debug_message("Removing model " . $self->__display_name__ . " from model group " . $group->__display_name__);
        my $rv = $group->unassign_models($self);
        unless ($rv) {
            Carp::confess $self->error_message("Failed to remove model " . $self->__display_name__ .
                " from model group " . $group->__display_name__);
        }
    }

    my @links = ($self->from_model_links, $self->to_model_links);
    for my $link (@links) {
        $self->debug_message("Deleting model link " . $link->__display_name__);
        my $rv = $link->delete;
        unless ($rv) {
            Carp::confess $self->error_message("Could not delete model link " . $link->__display_name__ .
                " prior to deleting model " . $self->__display_name__);
        }
    }

    # Get the remaining events like create and assign instrument data
    for my $event ($self->events) {
        unless ($event->delete) {
            $self->error_message('Failed to remove event '. $event->class .' '. $event->id);
            die $self->error_message();
        }
    }

    my $rv = $self->SUPER::delete(@_);
    return $rv;
}

sub set_apipe_cron_status {
    my $self = shift;
    my $body_text = shift;

    my @header = (header_text => 'apipe_cron_status');

    my $note = $self->notes(@header);
    $note->delete if ($note);

    $self->add_note(@header, body_text => $body_text);
}

sub duplicates {
    my $self    = shift || die;
    my $pp      = $self->processing_profile || die;
    my $class   = $self->class || die;
    my $subject = $self->subject || die;
    my @inputs  = $self->inputs;

    # duplicates would have the same subject, processing profile, and inputs
    # but we have to compare the values of the inputs not the inputs themselves
    my @duplicates;
    my @other_models;
    if (@_) {
        @other_models = grep { $_->subject_id eq $subject->id} @_;
    } else {
        @other_models = $class->get(subject_id => $subject->id, processing_profile_id => $pp->id);
    }

    my $instrument_data_ids = join(",", sort $self->instrument_data);
    for my $other_model (@other_models) {
        my $other_instrument_data_ids = join(",", sort $other_model->instrument_data);
        next unless $instrument_data_ids eq $other_instrument_data_ids;

        my @other_inputs = $other_model->inputs;
        next if (@other_inputs != @inputs); # mainly to catch case where one has inputs but other does not

        my $matched_inputs = 0;
        for my $input (@inputs) {
            my @other_duplicate_inputs = $other_model->inputs(name => $input->name, value_id => $input->value_id, value_class_name => $input->value_class_name);
            $matched_inputs++ if (@other_duplicate_inputs);
        }
        push @duplicates, $other_model if (@inputs == $matched_inputs);
    }

    @duplicates = grep { $_->id ne $self->id } @duplicates;

    return @duplicates;
}

sub is_used_as_model_or_build_input {
    # Both models and builds have this method and as such it is currently duplicated.
    # We don't seem to have any place to put things that are common between Models and Builds.
    my $self = shift;

    my @model_inputs = Genome::Model::Input->get(
        value_id => $self->id,
        value_class_name => $self->class,
    );

    my @build_inputs = Genome::Model::Build::Input->get(
        value_id => $self->id,
        value_class_name => $self->class,
    );

    my @inputs = (@model_inputs, @build_inputs);

    return (scalar @inputs) ? 1 : 0;
}

sub builds_are_used_as_model_or_build_input {
    my $self = shift;

    my @builds = $self->builds;

    return grep { $_->is_used_as_model_or_build_input } @builds;
}

sub dependent_properties {
    my ($self, $property_name) = @_;
    return;
}

sub _resolve_type_name_for_class {
    my $class = shift;
    my ($subclass) = $class =~ /^Genome::Model::([\w\d]+)$/;
    return unless $subclass;
    return Genome::Utility::Text::camel_case_to_string($subclass);
}

sub compatible_instrument_data {
    my $self = shift;
    my %params;

    $params{sequencing_platform} = $self->sequencing_platform if $self->sequencing_platform;
    if (my @samples = $self->get_all_possible_samples)  {
        my @sample_ids = map($_->id, @samples);
        %params = (
                   sample_id => \@sample_ids,
               );
    } else {
        %params = (
                   $self->subject_type => $self->subject_name,
               );
    }
    my @compatible_instrument_data = Genome::InstrumentData->get(%params);

    if($params{sequencing_platform} and $params{sequencing_platform} eq 'solexa') {
        # FASTQs with 0 reads crash in alignment.  Don't assign them. -??
        # TODO: move this into the assign logic, not here. -ss
        my @filtered_compatible_instrument_data;
        for my $idata (@compatible_instrument_data) {
            if (defined($idata->total_bases_read)&&($idata->total_bases_read == 0)) {
                $self->warning_message(sprintf("ignoring %s because it has zero bases read",$idata->__display_name__));
                next;
            }
            else {
                push @filtered_compatible_instrument_data, $idata;
            }
        }
        @compatible_instrument_data = @filtered_compatible_instrument_data;
    }

    return @compatible_instrument_data;
}

sub assigned_instrument_data { return $_[0]->instrument_data; }
sub available_instrument_data { return unassigned_instrument_data(@_); }
sub unassigned_instrument_data {
    my $self = shift;

    my @compatible_instrument_data = $self->compatible_instrument_data;
    my @assigned = $self->instrument_data;
    return @compatible_instrument_data unless @assigned;

    my %assigned_instrument_data_ids = map { $_->id => 1 } @assigned;
    return grep { not $assigned_instrument_data_ids{$_->id} } @compatible_instrument_data;
}
