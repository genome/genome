package Genome::Model::Command::Services::AssignQueuedInstrumentData;

use strict;
use warnings;

class Genome::Model::Command::Services::AssignQueuedInstrumentData {
    is  => 'Command::V2',
    has_optional => [
        instrument_data => {
            is          => 'Genome::InstrumentData',
            is_many     => 1,
            doc         => '[Re]process these instrument data.',
        },
        max_instrument_data_to_process => {
            is          => 'Number',
            default     => 200,
            doc         => 'Max # of instrument data to process in one invocation.',
        },
    ],
    has_transient => [
        _existing_models_with_existing_assignments => {
            is => 'HASH',
            doc => 'Existing models that already had the instrument data assigned.',
            default_value => {},
            is_output => 1,
        },
        _existing_models_assigned_to => {
            is => 'HASH',
            doc => 'Existing models with newly assigned instrument data.',
            default_value => {},
            is_output => 1,
        },
        _newly_created_models => {
            is => 'HASH',
            doc => 'New models created for the instrument data.',
            default_value => {},
            is_output => 1,
        },
    ],
};

sub _default_rna_seq_processing_profile_id {
    my $self = shift;
    my $instrument_data = shift;

    return 2762841;
}

sub help_brief {
    return 'Assign queued instrument data to models';
}

sub help_synopsis {
    return <<'EOS'
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;

    unless ( $ENV{UR_DBI_NO_COMMIT} ) {
        my $lock_resource = $ENV{GENOME_LOCK_DIR} . '/genome_model_command_services_assign-queued-instrument-data/loader';

        my $lock = Genome::Sys->lock_resource(resource_lock=>$lock_resource, max_try=>1);
        unless ($lock) {
            $self->error_message("could not lock, another instance must be running.");
            return;
        }

        UR::Context->current->add_observer(
            aspect => 'commit',
            callback => sub{
                Genome::Sys->unlock_resource(resource_lock=>$lock);
            }
        )
    }

    my @instrument_data_to_process = $self->_load_instrument_data;
    return 1 unless @instrument_data_to_process;

    INST_DATA: foreach my $instrument_data ( @instrument_data_to_process ) {
        $self->status_message('Starting instrument data '.$instrument_data->id);

        my $sequencing_platform = $instrument_data->sequencing_platform;

        my @processing = $self->_resolve_processing_for_instrument_data($instrument_data);

        if ( $instrument_data->ignored ) { # could be set in _resolve_processing_for_instrument_data
            $self->status_message('Skipping ignored instrument data! '.$instrument_data->id);
            next INST_DATA;
        }

        my $subject = $instrument_data->sample;
        my @process_errors;
        if ( @processing ) {
            PP: foreach my $processing ( @processing ) {
                if(exists $processing->{error}) {
                    push @process_errors, $processing->{error};
                    next PP;
                }

                my $processing_profile = $processing->{processing_profile};
                my $reference_sequence_build = $processing->{reference_sequence_build};

                # These pps require reference seq build
                my @processing_profiles_that_require_reference_sequence_build = (qw/
                    Genome::ProcessingProfile::ReferenceAlignment
                    Genome::ProcessingProfile::GenotypeMicroarray
                    Genome::ProcessingProfile::RnaSeq
                    /);
                if ( not $reference_sequence_build and grep { $processing_profile->isa($_) } @processing_profiles_that_require_reference_sequence_build ) {
                    $self->error_message('Failed to set reference sequence build during processing instrument data! '.$instrument_data->id);
                    push @process_errors, $self->error_message;
                    next PP;
                }

                my @models = Genome::Model->get(
                    subject_id            => $subject->id,
                    processing_profile_id => $processing_profile->id,
                    auto_assign_inst_data => 1,
                );

                my @assigned = $self->assign_instrument_data_to_models($instrument_data, $reference_sequence_build, @models);

                #returns an explicit undef on error
                if(scalar(@assigned) eq 1 and not defined $assigned[0]) {
                    push @process_errors, $self->error_message;
                    next PP;
                }


                if(scalar(@assigned > 0)) {
                    #find or create default qc models if applicable
                    $self->create_default_qc_models(@assigned);
                    #find or create somatic models if applicable
                    $self->find_or_create_somatic_variation_models(@assigned);

                } else {
                    # no model found for this PP, make one (or more) and assign all applicable data
                    my @new_models = $self->create_default_models_and_assign_all_applicable_instrument_data($instrument_data, $subject, $processing_profile, $reference_sequence_build);
                    unless(@new_models) {
                        push @process_errors, $self->error_message;
                        next PP;
                    }
                    #find or create somatic models if applicable
                    $self->find_or_create_somatic_variation_models(@new_models);
                }
            } # looping through processing profiles for this instdata, finding or creating the default model
        } elsif ( $sequencing_platform eq 'solexa'
                and $instrument_data->target_region_set_name
                and Genome::FeatureList->get(name => $instrument_data->target_region_set_name)
                and Genome::FeatureList->get(name => $instrument_data->target_region_set_name)->content_type eq 'validation'
        ) {
            my @validation = Genome::Model::SomaticValidation->get(
                target_region_set_name => $instrument_data->target_region_set_name,
            );

            @validation = grep((($_->tumor_sample and $_->tumor_sample eq $instrument_data->sample) or ($_->normal_sample and $_->normal_sample eq $instrument_data->sample)), @validation);
            if(@validation) {
                my $fl = Genome::FeatureList->get(name => $instrument_data->target_region_set_name);
                my $ok = 0;

                #try all possible matching references
                for($fl->reference, map($_->destination_reference_build, Genome::Model::Build::ReferenceSequence::Converter->get(source_reference_build_id => $fl->reference->id)) ) {
                    $ok = $self->assign_instrument_data_to_models($instrument_data, $_, @validation) || $ok;
                }

                unless($ok) {
                    push @process_errors,
                    $self->error_message('Did not assign validation instrument data to any models.');
                }
            } elsif($instrument_data->index_sequence eq 'unknown' && $instrument_data->sample->name =~ /Pooled_Library/) {
                $self->status_message('Skipping pooled library validation data! '.$instrument_data->id);
                $instrument_data->{_skipped} = 1;
            } else {
                push @process_errors,
                $self->error_message('No validation models found to assign data (target ' . $instrument_data->target_region_set_name . ' on instrument data ' . $instrument_data->id . '.)');
            }
        } else {
            $self->status_message('No model generation attempted for instrument data! '.$instrument_data->id);
            $instrument_data->{_skipped} = 1;
        } # done with inst data which specify @processing

        # Handle this instdata for other models besides the default
        {
            my @found_models;
            my @check = qw/sample taxon/;

            for my $check (@check) {
                my $subject = $instrument_data->$check;
                if (defined($subject)) {
                    my @some_models= Genome::Model->get(
                        subject_id         => $subject->id,
                        auto_assign_inst_data => 1,
                    );

                    my $new_models = $self->_newly_created_models;
                    @some_models = grep { not $new_models->{$_->id} } @some_models;
                    push @found_models,@some_models;
                }
            }

            @found_models =
            grep {
                $_->processing_profile->can('sequencing_platform')
            } @found_models;

            @found_models =
            grep {
                $_->processing_profile->sequencing_platform() eq $sequencing_platform
            } @found_models;

            #Don't care here what ref. seq. was used (if any)
            my @assigned = $self->assign_instrument_data_to_models($instrument_data, undef, @found_models);
            if(scalar(@assigned) eq 1 and not defined $assigned[0]) {
                push @process_errors, $self->error_message;
            }
        } # end of adding instdata to non-autogen models

        $instrument_data->{_tgi_lims_fail_message} = substr(join("\n", @process_errors), 0, 512) if @process_errors; # max for msg is 512
    } # end of INST_DATA loop

    #schedule new builds for the models we found and stored in the output hashes
    $self->request_builds;

    $self->status_message("Updating instrument data...");
    for my $instrument_data ( @instrument_data_to_process ) {
        if ( $instrument_data->ignored or $instrument_data->{_skipped} ) { # skipped
            $self->_update_instrument_data_tgi_lims_status_to_skipped($instrument_data);
        }
        elsif ( $instrument_data->{_tgi_lims_fail_message} ) { # failed
            $self->_update_instrument_data_tgi_lims_status_to_failed($instrument_data);
        }
        else { # processed!
            $self->_update_instrument_data_tgi_lims_status_to_processed($instrument_data);
        }
    }

    return 1;
}

sub _update_instrument_data_tgi_lims_status_to {
    my ($self, $instrument_data, $status) = @_;

    # These should not happen - developer error
    Carp::confess('No instrument data given to update instrument data tgi lims status!') if not $instrument_data;
    Carp::confess('No status given to update instrument data tgi lims status!') if not $status;
    Carp::confess("No invalid status ($status) given to update instrument data tgi lims status!") if not grep { $status eq $_ } (qw/ processed skipped failed /);

    # Rm tgi lims status attribute(s)
    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_status');

    # Set tgi lims status attribute to processed
    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_status',
        attribute_value => $status,
    );
    $self->status_message("Set TGI LIMS status to $status for instrument data! ".$instrument_data->id);

    return 1;
}

sub _update_instrument_data_tgi_lims_status_to_processed {
    my ($self, $instrument_data) = @_;

    my $set_status = $self->_update_instrument_data_tgi_lims_status_to($instrument_data, 'processed');
    return if not $set_status;

    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_message');
    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_count');

    return 1;
}

sub _update_instrument_data_tgi_lims_status_to_skipped {
    my ($self, $instrument_data) = @_;

    my $set_status = $self->_update_instrument_data_tgi_lims_status_to($instrument_data, 'skipped');
    return if not $set_status;

    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_message');
    $instrument_data->remove_attribute(attribute_label => 'tgi_lims_fail_count');

    return 1;
}

sub _update_instrument_data_tgi_lims_status_to_failed {
    my ($self, $instrument_data) = @_;

    Carp::confess('No message sent to update tgi lims status to failed!') if not defined $instrument_data->{_tgi_lims_fail_message};

    my $set_status = $self->_update_instrument_data_tgi_lims_status_to($instrument_data, 'failed');
    return if not $set_status;

    my $fail_msg_attr = $instrument_data->attributes(attribute_label => 'tgi_lims_fail_message');
    $fail_msg_attr->delete if $fail_msg_attr;

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_fail_message',
        attribute_value => $instrument_data->{_tgi_lims_fail_message},
    );

    my $fail_count_attr = $instrument_data->attributes(attribute_label => 'tgi_lims_fail_count');
    my $previous_count = 0;
    if ( $fail_count_attr ) {
        $previous_count = $fail_count_attr->attribute_value;
        $fail_count_attr->delete;
    }

    $instrument_data->add_attribute(
        attribute_label => 'tgi_lims_fail_count',
        attribute_value => ($previous_count+1),
    );

    return 1;
}

sub find_or_create_somatic_variation_models{
    my ($self, @models) = @_;
    #only want sample-based models
    @models = grep { $_->subject_class_name eq 'Genome::Sample' } @models;
    #only want TCGA models
    @models = grep {$self->is_tcga_reference_alignment($_) } @models;
    #We want capture models with one of the given roi_set_names and all non capture models here.  Filter the rest out
    @models = grep {defined($_->region_of_interest_set_name) ? $_->region_of_interest_set_name =~ m/agilent.sureselect.exome.version.2.broad.refseq.cds.only/ : 1} @models;
    for my $model (@models){
        my $sample = $model->subject;
        #find or create mate ref-align model
        if ($sample->name =~ m/([^-]+-[^-]+-[^-]+-)([01]{2}A)(.*)/){
            my ($prefix, $designator, $suffix) = ($1, $2, $3);
            my %designator_pairing = (
                '10A' => '01A',
                '01A' => '10A',
            );
            my $mate_designator = $designator_pairing{$designator};
            $self->error_message("Not processing somatic variation model with sample name: " . $sample->name . " and designator: $designator") and next unless $mate_designator;
            my $mate_name = join("", $prefix, $mate_designator, $suffix);

            my $subject_for_mate = Genome::Sample->get(name => $mate_name);
            $self->error_message("No sample found for mate_name $mate_name (paired to: " . $model->name . ")") and next unless $subject_for_mate;

            my %mate_params = (
                subject_id => $subject_for_mate->id,
                reference_sequence_build => $model->reference_sequence_build,
                processing_profile => $model->processing_profile,
                auto_assign_inst_data => '1',
            );
            $mate_params{annotation_reference_build_id} = $model->annotation_reference_build_id if $model->can('annotation_reference_build_id') and $model->annotation_reference_build_id;
            $mate_params{target_region_set_name} = $model->target_region_set_name if $model->can('target_region_set_name') and $model->target_region_set_name;
            $mate_params{region_of_interest_set_name} = $model->region_of_interest_set_name if $model->can('region_of_interest_set_name') and $model->region_of_interest_set_name;

            my ($mate) = Genome::Model::ReferenceAlignment->get( %mate_params );
            unless ($mate){
                $mate = $model->copy(
                    name => 'AQID-PLACE_HOLDER',
                    instrument_data => undef,
                );
                $self->error_message("Failed to find copied mate with subject name: $mate_name") and next unless $mate;

                $mate->subject_id($subject_for_mate->id);

                my $capture_target = eval{$model->target_region_set_name};
                my $mate_model_name = $mate->default_model_name(capture_target => $capture_target);
                $self->error_message("Could not name mate model for with subject name: $mate_name") and next unless $mate_model_name;
                $mate->name($mate_model_name);
                $mate->auto_assign_inst_data(1);
                $mate->build_requested(0, 'AQID: newly created mate for creating somatic-variation model--has no instrument data');

                my $new_models = $self->_newly_created_models;
                $new_models->{$mate->id} = $mate;
            }

            my %somatic_params = (
                auto_assign_inst_data => 1,
            );
            $somatic_params{annotation_build} = Genome::Model::ImportedAnnotation->annotation_build_for_reference($model->reference_sequence_build);
            $self->error_message('Failed to get annotation_build for somatic variation model with model: ' . $model->name) and next unless $somatic_params{annotation_build};
            $somatic_params{previously_discovered_variations_build} = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($model->reference_sequence_build);
            $self->error_message('Failed to get previously_discovered_variations_build for somatic variation model with model: ' . $model->name) and next unless $somatic_params{previously_discovered_variations_build};
 

            my $capture_somatic_processing_profile_id = '2762563'; #Oct 2012 Default Somatic Variation Exome
            my $somatic_processing_profile_id = '2762562'; #Oct 2012 Default Somatic Variation WGS
            my $capture_target = eval{$model->target_region_set_name};
            if($capture_target){
                $somatic_params{processing_profile_id} = $capture_somatic_processing_profile_id;
            }
            else{
                $somatic_params{processing_profile_id} = $somatic_processing_profile_id;
            }
            if ($designator eq '10A'){ #$model is normal
                $somatic_params{normal_model} = $model;
                $somatic_params{tumor_model} = $mate;
            }elsif ($designator eq '01A'){ #$model is tumor
                $somatic_params{tumor_model} = $model;
                $somatic_params{normal_model} = $mate;
            }else{
                die $self->error_message("Serious error in sample designators for automated create of somatic-variation models for ".$model->subject_name);
            }

            my $somatic_variation = Genome::Model::SomaticVariation->get(%somatic_params);

            unless ($somatic_variation){
                $somatic_params{model_name} = 'AQID-PLACE_HOLDER';
                my $create = Genome::Model::Command::Define::SomaticVariation->execute( %somatic_params );
                $self->error_message('Failed to create somatic variation model with component model: ' . $model->name) and next unless $create;

                delete $somatic_params{model_name};
                $somatic_params{name} = 'AQID-PLACE_HOLDER';
                $somatic_variation = Genome::Model::SomaticVariation->get(%somatic_params);
                $self->error_message("Failed to find new somatic variation model with component model: " . $model->name) and next unless $somatic_variation;

                $somatic_variation->build_requested(0, 'AQID: somatic variation build is not ready until ref. align. builds finish');
                my $somatic_variation_model_name = $somatic_variation->default_model_name(capture_target => $capture_target);
                $self->error_message("Failed to name new somatic variation model with component model: " . $model->name) and next unless $somatic_variation_model_name;
                $somatic_variation->name($somatic_variation_model_name);

                my $new_models = $self->_newly_created_models;
                $new_models->{$somatic_variation->id} = $somatic_variation;
            }
        }
        else{
            $self->error_message("Not processing somatic variation model with sample name: " . $model->subject_name);
        }


    }
}

sub first_compatible_feature_list_name {
    my $self = shift || die;
    my $reference = shift || die;
    my $feature_list_names = shift || die;
    for my $feature_list_name (@$feature_list_names) {
        my $feature_list = Genome::FeatureList->get(name => $feature_list_name);
        unless ($feature_list) {
            croak("non-existant FeatureList name ($feature_list_name) passed");
        }

        my $feature_list_reference = $feature_list->reference;
        unless ($feature_list_reference) {
            croak("reference not set for FeatureList (name = $feature_list_name)");
        }

        if ($reference->is_compatible_with($feature_list_reference)) {
            return $feature_list_name;
        }
    }

    return;
}

sub is_tcga_reference_alignment {
    my $self = shift;
    my $model = shift;
    my $sample = $model->subject;

    return unless $model->isa('Genome::Model::ReferenceAlignment');
    return if ($model->isa('Genome::Model::ReferenceAlignment') && $model->is_lane_qc);

    #try the extraction label
    my @results = grep {$_->attribute_label eq 'extraction_label' and $_->attribute_value =~ m/^TCGA/} $sample->attributes;
    return 1 if @results;

    #otherwise, check the nomenclature
    my @nomenclature = map { $_->nomenclature } ($sample, $sample->attributes);
    return grep { $_ && $_ =~ /^TCGA/i } @nomenclature;
}

sub _load_instrument_data {
    my $self = shift;

    $self->status_message('Get instrument data...');
    
    my %instrument_data;
    if ( $self->instrument_data ) {
        %instrument_data = map { $_->{_priority} = 0; $_->id => $_ } $self->instrument_data;
    }
    else {
        my @status_attrs = Genome::InstrumentDataAttribute->get(
            attribute_label => 'tgi_lims_status',
            attribute_value => [qw/ new failed /],
        );
        for my $status_attr ( @status_attrs ) {
            my $instrument_data = Genome::InstrumentData->get(
                id => $status_attr->instrument_data_id,
                analysis_project_id => undef,
                -hint => [ 'sample', 'sample.source', 'sample.source.taxon', ],
            );
            if ($instrument_data) {
                my $fail_cnt = eval{ $instrument_data->attributes(attribute_label => 'tgi_lims_fail_count')->attribute_value; };
                $instrument_data->{_priority} = ( $fail_cnt ? $fail_cnt : 0 ); # if it does not have a fail count, treat as new
                $instrument_data{ $instrument_data->id } = $instrument_data;
            }
        }
    }
    $self->status_message('Found '.scalar(grep { $_->{_priority} == 0 } values %instrument_data)." new instrument data\n");
    $self->status_message('Found '.scalar(grep { $_->{_priority} > 0 } values %instrument_data)." previously attempted instrument data\n");

    $self->status_message('Filter instrument data we can process...');
    # sort for determinism since data is in a hash
    my $sorter = sub { $a->{_priority} <=> $b->{_priority} or Genome::InstrumentData->__meta__->id_property_sorter->($a, $b) };
    my @instrument_data_to_process;
    my $max_instrument_data_to_process = $self->max_instrument_data_to_process;
    for my $instrument_data ( sort $sorter values %instrument_data ) {
        last if @instrument_data_to_process >= $max_instrument_data_to_process;
        if ( not $self->_check_instrument_data($instrument_data) ){
            $self->_update_instrument_data_tgi_lims_status_to_failed($instrument_data);
            next;
        }
        push @instrument_data_to_process, $instrument_data;
    }
    $self->status_message('Processing '.@instrument_data_to_process.' instrument data');

    return @instrument_data_to_process;
}

sub _check_instrument_data {
    my ($self, $instrument_data) = @_;

    if ( $instrument_data->isa('Genome::InstrumentData::Solexa') ) {
        if($instrument_data->target_region_set_name) {
            my $fl = Genome::FeatureList->get(name => $instrument_data->target_region_set_name);
            if(not $fl) {
                $instrument_data->{_tgi_lims_fail_message} = 'Failed to get a feature-list matching target region set name ' . $instrument_data->target_region_set_name;
            }
            elsif(not $fl->content_type) {
                $instrument_data->{_tgi_lims_fail_message} = 'No content-type set on feature-list ' . $fl->name;
            } elsif ($fl->content_type eq 'roi') {
                $instrument_data->{_tgi_lims_fail_message} = 'Unexpected "roi"-typed feature-list set as target region set name: ' . $fl->name;
            } elsif (!grep($_ eq $fl->content_type, 'exome', 'validation', 'targeted')) {
                $instrument_data->{_tgi_lims_fail_message} = 'Unknown/unhandled content-type ' . $fl->content_type . ' on feature-list ' . $fl->name;
            }
        }
    }

    if ( $instrument_data->{_tgi_lims_fail_message} ) {
        $self->error_message($instrument_data->{_tgi_lims_fail_message});
        return;
    }

    return 1;
}

sub assign_instrument_data_to_models {
    my $self = shift;
    my $genome_instrument_data = shift;
    my $reference_sequence_build = shift;
    my @models = @_;

    my $instrument_data_id = $genome_instrument_data->id;

    #only assign to models that have auto_assign_inst_data=1
    @models = grep($_->auto_assign_inst_data, @models);

    # we don't want to (automagically) assign capture and non-capture data to the same model.
    if ( @models and $genome_instrument_data->can('target_region_set_name') ) {
        my $id_capture_target = $genome_instrument_data->target_region_set_name();

        if ($id_capture_target) {
            # keep only models with the specified capture target
            @models = grep($_->can('target_region_set_name') && $_->target_region_set_name && $_->target_region_set_name eq $id_capture_target, @models);
        } else {
            # keep only models with NO capture target
            my %capture_model_ids = map { $_->model_id => 1 } Genome::Model::Input->get(
                model_id => [ map { $_->id } @models ],
                name => 'target_region_set_name',
            );
            @models = grep { not $capture_model_ids{$_->id} } @models;
        }
    }

    #we don't want to (automagically) assign rna seq and non-rna seq data to the same model.
    unless ($genome_instrument_data->isa('Genome::InstrumentData::454')) { #454 data should be allowed to be in MC16S models that it's explicitly looking for
        if (@models and $genome_instrument_data->sample->is_rna) {
            @models = grep($_->isa('Genome::Model::RnaSeq'), @models);
        }else{
            @models = grep(!($_->isa('Genome::Model::RnaSeq')), @models);
        }
    }

    if($reference_sequence_build) {
        @models = grep($_->reference_sequence_build eq $reference_sequence_build, @models);
    }

    foreach my $model (@models) {
        my @existing_instrument_data = $model->input_for_instrument_data_id($instrument_data_id);

        if (@existing_instrument_data) {
            $self->warning_message(
                "instrument data '$instrument_data_id'" .
                ' already assigned to model ' . $model->id
            );

            my $existing_models = $self->_existing_models_with_existing_assignments;
            $existing_models->{$model->id} = $model;
        } else {
            my $assign =
            Genome::Model::Command::InstrumentData::Assign->create(
                instrument_data => [$genome_instrument_data],
                model           => $model,
            );

            unless ( $assign->execute ) {
                $self->error_message(
                    'Failed to execute instrument-data assign for '
                    . 'model '
                    . $model->id
                    . ' and instrument data '
                    . $instrument_data_id );
                return undef;
            }

            my $existing_models = $self->_existing_models_assigned_to;
            $existing_models->{$model->id} = $model;
        }
    }
    return @models;
}

sub create_default_models_and_assign_all_applicable_instrument_data {
    my $self = shift;
    my $genome_instrument_data = shift;
    my $subject = shift;
    my $processing_profile = shift;
    my $reference_sequence_build = shift;

    my @new_models;
    my @ref_align_models;

    my %model_params = (
        name                    => 'AQID-PLACE_HOLDER',
        user_name               => 'apipe-builder',
        subject_id              => $subject->id,
        subject_class_name      => $subject->class,
        processing_profile_id   => $processing_profile->id,
        auto_assign_inst_data   => 1,
    );

    if ( grep { my $class = 'Genome::ProcessingProfile::'.$_; $processing_profile->isa($class); } (qw/ GenotypeMicroarray RnaSeq /) ) {
        $model_params{auto_assign_inst_data} = 0;
    }

    if ( $reference_sequence_build ) {
        $model_params{reference_sequence_build} = $reference_sequence_build;
        unless( $processing_profile->isa('Genome::ProcessingProfile::RnaSeq')){
            my $dbsnp_build = Genome::Model::ImportedVariationList->dbsnp_build_for_reference($reference_sequence_build);
            $model_params{dbsnp_build} = $dbsnp_build if $dbsnp_build;
        }

        #annotion build inputs
        if ( $processing_profile->isa('Genome::ProcessingProfile::ReferenceAlignment')){
            my $annotation_build = Genome::Model::ImportedAnnotation->annotation_build_for_reference($reference_sequence_build);
            $model_params{annotation_reference_build} = $annotation_build if $annotation_build;
        }
        if ( $processing_profile->isa('Genome::ProcessingProfile::RnaSeq')){
            my $annotation_build = Genome::Model::ImportedAnnotation->annotation_build_for_reference($reference_sequence_build);
            $model_params{annotation_build} = $annotation_build if $annotation_build;
        }
    }

    my $regular_model = Genome::Model->create(%model_params);
    unless ( $regular_model ) {
        $self->error_message('Failed to create model with params: '.Data::Dumper::Dumper(\%model_params));
        return;
    }
    push @new_models, $regular_model;

    my $capture_target = eval{ $genome_instrument_data->target_region_set_name; };

    my $name = $regular_model->default_model_name(
        instrument_data => $genome_instrument_data,
        capture_target => $capture_target,
    );
    if ( not $name ) {
        $self->error_message('Failed to get model name for params: '.Data::Dumper::Dumper(\%model_params));
        for my $model ( @new_models ) { $model->delete; }
        return;
    }
    $regular_model->name($name);

    if ($regular_model->isa('Genome::Model::ReferenceAlignment')) {
        push @ref_align_models, $regular_model;
    }

    if ( $capture_target and $reference_sequence_build and not $regular_model->isa('Genome::Model::RnaSeq')){
        # FIXME This is a lame hack for these capture sets
        my %build36_to_37_rois = get_build36_to_37_rois();
        my $root_build37_ref_seq = Genome::Model::Build::ImportedReferenceSequence->get(106942997) || die;# GRCh37-lite-build37 => 106942997
        if ( not $root_build37_ref_seq ) {
            $self->error_message('Failed to get reference sequence build for id! '.106942997);
            return;
        }

        my $roi_list = $capture_target;
        if ($reference_sequence_build
                and $reference_sequence_build->is_compatible_with($root_build37_ref_seq)
                and exists $build36_to_37_rois{$capture_target}
        ) {
            $roi_list = $build36_to_37_rois{$capture_target};
        }

        unless($self->assign_capture_inputs($regular_model, $capture_target, $roi_list)) {
            for my $model ( @new_models ) { $model->delete; }
            return;
        }

        my %roi_sets = roi_sets();

        for my $roi_set (keys %roi_sets) {
            my $roi_set_names = $roi_sets{$roi_set};
            my $roi_set_name = $self->first_compatible_feature_list_name($reference_sequence_build, $roi_set_names);
            if ($roi_set_name) {
                my $roi_set_model = $self->create_roi_model($roi_set, $genome_instrument_data, $roi_set_name, %model_params);
                if ($roi_set_model) {
                    push @new_models, $roi_set_model;
                } else {
                    $self->error_message("Failed to create $roi_set model.");
                    for my $model (@new_models) { $model->delete; }
                    return;
                }
            }
        }
    }

    for my $m (@new_models) {
        my $assign =
        Genome::Model::Command::InstrumentData::Assign->create(
            model => $m,
            instrument_data => [$genome_instrument_data],
            include_imported => 1,
            force => 1,
        );

        unless ( $assign->execute ) {
            $self->error_message(
                'Failed to execute instrument-data assign for model '
                . $m->id . ' instrument data '.$genome_instrument_data->id );
            for my $model (@new_models) { $model->delete; }
            return;
        }

        unless($m->isa('Genome::Model::RnaSeq')){
            my $assign_all =
            Genome::Model::Command::InstrumentData::Assign->create(
                model => $m,
                all => 1,
            );

            unless ( $assign_all->execute ) {
                $self->error_message(
                    'Failed to execute instrument-data assign --all for model '
                    . $m->id );
                for my $model (@new_models) { $model->delete; }
                return;
            }
        }

        my @existing_instrument_data = $m->input_for_instrument_data($genome_instrument_data);
        unless (@existing_instrument_data) {
            $self->error_message(
                'instrument data ' . $genome_instrument_data->id . ' not assigned to model ????? (' . $m->id . ')'
            );
            for my $model (@new_models) { $model->delete; }
            return;
        }

        unless($m->isa('Genome::Model::GenotypeMicroarray')){
            $self->add_model_to_default_projects($m, $genome_instrument_data);
        }

        my $new_models = $self->_newly_created_models;
        $new_models->{$m->id} = $m;
    }

    # Now that they've had their instrument data assigned get_or_create_lane_qc_models
    # Based of the ref-align models so that alignment can shortcut
    push(@new_models , $self->create_default_qc_models(@ref_align_models));
    return @new_models;
}

sub create_roi_model {
    my $self = shift;
    my $roi = shift || die;
    my $genome_instrument_data = shift || die;
    my $region_of_interest_set_name = shift || die;
    my %model_params = @_;
    die unless keys %model_params;

    my $abortion_message = sub {
        my $message = shift;
        return $self->error_message("Aborting creation of $roi model, $message.");
    };
    my $deletion_message = sub {
        my $message = shift;
        return $self->error_message("Deleting partially created $roi model, $message.");
    };

    my $target_region_set_name = eval { $genome_instrument_data->target_region_set_name };
    unless ($target_region_set_name) {
        $abortion_message->('instrument data does not have a target region');
        return;
    }

    my $model = Genome::Model->create(%model_params);
    unless ($model) {
        $abortion_message->('failed to create model');
        return;
    }

    my $roi_model_name = $model->default_model_name(
        instrument_data => $genome_instrument_data,
        capture_target => $target_region_set_name,
        roi => lc($roi),
    );
    unless ($roi_model_name) {
        $deletion_message->('failed to resolve default model name');
        $model->delete;
        return;
    }

    unless ($model->name($roi_model_name)) {
        $deletion_message->('failed to rename to default model name');
        $model->delete;
        return;
    }

    unless($self->assign_capture_inputs($model, $target_region_set_name, $region_of_interest_set_name)) {
        $deletion_message->('failed to assign capture inputs');
        $model->delete;
        return;
    }

    return $model;
}

sub get_build36_to_37_rois {
    return (
        'agilent sureselect exome version 2 broad refseq cds only' => 'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37',
        'agilent sureselect exome version 2 broad' => 'agilent sureselect exome version 2 broad hg19 liftover',
        'hg18 nimblegen exome version 2' => 'hg19 nimblegen exome version 2',
        'NCBI-human.combined-annotation-54_36p_v2_CDSome_w_RNA' => 'NCBI-human.combined-annotation-54_36p_v2_CDSome_w_RNA_build36-build37_liftOver',
        'Freimer Pool of original (4k001L) plus gapfill (4k0026)' => 'Freimer-Boehnke capture-targets.set1_build37-fix1',
        '04110401 PoP32 EZ capture chip set'   => '04110401 PoP32 EZ capture chip set build37',
        'RT 49315 - AMD -- pool 1' => 'AMD-pool1-build37',
        '03110401 capture chip set' => '03110401 capture chip set - liftover_build37',
        'CleftPalate 03110402 capture chip set' => 'CleftPalate 03110402 capture chip set - liftover_build37',
        '34010 capture oligo tube' => '34010 capture oligo tube -- liftover_build37',
    );
}

sub create_default_qc_models {
    my $self = shift;
    my @models = @_;
    my @new_models;
    for my $model (@models){
        next unless $model->type_name eq 'reference alignment';
        next unless $model->processing_profile_name =~ /^\w+\ \d+\ Default\ Reference\ Alignment/; # e.g. Feb 2011 Defaulte Reference Alignment

        my @lane_qc_models = $model->get_or_create_lane_qc_models;

        for my $lane_qc (@lane_qc_models) {
            next if $lane_qc->build_requested;
            next unless $lane_qc->build_needed;
            $lane_qc->build_requested(1);
            push @new_models, $lane_qc;
        }
    }

    return @new_models;
}

sub _find_or_create_mc16s_454_qc_model {
    my ($self, $instrument_data) = @_;

    $self->status_message("Find or create mc16s 454 qc model!");

    my $pp_id = Genome::Model::MetagenomicComposition16s->default_processing_profile_id;
    my $name = $instrument_data->run_name.'_r'.$instrument_data->region_number.'.prod-mc16s-qc';
    my $model = Genome::Model->get(
        name => $name,
        processing_profile_id => $pp_id,
    );
    if ( not $model ) {
        my $model = Genome::Model->create(
            name => $name,
            subject_id => 2863615589, # Human Metagenome
            subject_class_name => 'Genome::Taxon',
            processing_profile_id => $pp_id,
            auto_assign_inst_data => 0,
        );
        $model->add_instrument_data($instrument_data);
        my $new_models = $self->_newly_created_models;
        $new_models->{$model->id} = $model;
    }
    else {
        my $existing_instrument_data = $model->inputs(name => 'instrument_data', value => $instrument_data);
        if ( not $existing_instrument_data ) {
            $model->add_instrument_data($instrument_data);
            $self->_existing_models_assigned_to->{$model->id} = $model;
        }
        else {
            $self->_existing_models_with_existing_assignments->{$model->id} = $model;
        }
    }

    return 1;
}

sub assign_capture_inputs {
    my $self = shift;
    my $model = shift;
    my $target_region_set_name = shift;
    my $region_of_interest_set_name = shift;

    my $target_input = $model->add_input(
        name             => "target_region_set_name",
        value_class_name => "UR::Value",
        value_id         => $target_region_set_name
    );

    unless ( defined($target_input) ) {
        $self->error_message('Failed to set capture target input for model ' . $model->id);
        return;
    }

    my $roi_input = $model->add_input(
        name             => "region_of_interest_set_name",
        value_class_name => "UR::Value",
        value_id         => $region_of_interest_set_name
    );

    my $fl = Genome::FeatureList->get(name => $region_of_interest_set_name);
    if($fl->content_type eq 'targeted' and $fl->is_multitracked and $model->isa('Genome::Model::ReferenceAlignment')) {
        $model->add_input(name => 'roi_track_name', value_class_name => 'UR::Value', value_id => 'tiled_region');
    }

    unless (defined($roi_input)) {
        $self->error_message('Failed to set region of instrument input for model ' . $model->id);
        return;
    }

    return 1;
}

sub _get_name_ext_for_model {
    my $self = shift;
    my $model = shift;

    my %roi_sets = roi_sets();
    for my $ext (keys %roi_sets) {
        my @roi_names = @{$roi_sets{$ext}};
        for my $roi_name (@roi_names) {
            return lc('.'.$ext) if $model->region_of_interest_set_name =~ /$roi_name/i;
        }
    }
}

sub add_model_to_default_projects {
    my $self = shift;
    my $model = shift;
    my $instrument_data = shift;

    # Get projects associated with the instrument data
    my @projects = $self->_get_projects_for_instrument_data($instrument_data);
    my $ext = $self->_get_name_ext_for_model($model);

    if ( $ext ) { # Get the projects for these names, but with the ext
        @projects = map { $self->_get_or_create_project_by_name($_) } map { $_->name.$ext } @projects;
    }
    for my $project ( @projects ) {
        $project->add_part(entity => $model);
    }

    # Get/create source project
    my $source_project_name = eval { # eval to return the name/undef
        my $subject = $model->subject;
        my $source;
        if ( $subject->isa('Genome::Sample') ) {
            $source = $subject->source;
        } elsif ( $subject->isa('Genome::Individual') ) {
            $source = $subject;
        } elsif ( $subject->isa('Genome::Library')) {
            my $sample = $subject->sample;
            $source = $sample->source;
        }
        return if not $source;

        my $common_name = $source->common_name;
        return if not $common_name;
        $common_name =~ /^([a-z]+)\d+$/i;
        return $1;
    };
    if ( $source_project_name ) {
        my $source_project = $self->_get_or_create_project_by_name($source_project_name);# dies on fail
        $source_project->add_part(entity => $model);
    }

    # Get/create pooled sample projects
    if ( my $pooled_sample_name = $self->_resolve_pooled_sample_name_for_instrument_data($instrument_data) ) {
        if ($ext =~ /\.wu-space$/) {
            $pooled_sample_name .= ".wu-space";
        }
        if ($ext =~ /\.tcga-cds$/) {
            $pooled_sample_name .= ".tcga-cds";
        }
        my $pooled_sample_project = $self->_get_or_create_project_by_name($pooled_sample_name); # dies on fail
        $pooled_sample_project->add_part(entity => $model);
    }

    return 1;
}

sub _get_projects_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my @project_parts = Genome::ProjectPart->get(
        entity_id => $instrument_data->id,
        label => 'instrument_data',
    );
    return if not @project_parts;

    my @projects = Genome::Project->get(
        id => [ map { $_->project_id } @project_parts ],
    );
    if ( not @projects or @projects != @project_parts ) {
        die $self->error_message('Failed to get projects for parts! '.Data::Dumper::Dumper(\@project_parts));
    }

    return @projects
}

sub _get_or_create_project_by_name {
    my ($self, $name) = @_;

    my $project = Genome::Project->get(name => $name);
    if ( not $project ) {
        $project = Genome::Project->create(
            name => $name,
        );
        if ( not $project ) {
            die $self->error_message('Failed to create project for name! '.$name);
        }
    }

    return $project;
}

sub _get_default_processing_profile_ids_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my @projects = $self->_get_projects_for_instrument_data($instrument_data);
    return if not @projects;

    my %pps;
    for my $project ( @projects ) {
        my @parts = $project->parts( label => 'default_processing_profiles' );
        next if not @parts;
        for my $part ( @parts ) {
            $pps{ $part->entity_id } = 1;
        }
    }

    return keys %pps;
}

sub _resolve_pooled_sample_name_for_instrument_data {
    my $self = shift;
    my $instrument_data = shift;

    return unless $instrument_data->can('index_sequence');
    my $index = $instrument_data->index_sequence;
    if($index) {
        my $instrument_data_class = $instrument_data->class;
        my $pooled_subset_name = $instrument_data->subset_name;
        $pooled_subset_name =~ s/${index}$/unknown/;

        my $pooled_instrument_data = $instrument_data_class->get(
            run_name => $instrument_data->run_name,
            subset_name => $pooled_subset_name,
            index_sequence => 'unknown',
        );
        return unless $pooled_instrument_data;

        my $sample = $pooled_instrument_data->sample;
        return unless $sample;

        return $sample->name;
    }

    return;
}

sub request_builds {
    my $self = shift;

    my $new_models = $self->_newly_created_models;
    my $assigned_to = $self->_existing_models_assigned_to;
    my %models_to_build;
    for my $model (values %$new_models) {
        #some models are explicitly not being built right away
        #but they might be picked up in other categories if instrument data is picked up in same AQID run
        next if defined $model->build_requested;
        $models_to_build{$model->id} = [$model, 'it is newly created'];
    }
    for my $model (values %$assigned_to) {
        next if exists $models_to_build{$model->id}; #already added above
        $models_to_build{$model->id} = [$model, 'it has been assigned to'];
    }

    $self->status_message("Finding models which need to build...");
    my $possibly_build = ($self->_existing_models_with_existing_assignments);
    for my $model (values %$possibly_build) {
        next if exists $models_to_build{$model->id}; #already added above

        my $last_build = $model->latest_build;

        unless(defined $last_build) {
            #no builds--can't possibly have built with all data
            my $reason = 'it has no builds';
            $self->status_message('Requesting build of model ' . $model->__display_name__ . " because $reason.");
            $models_to_build{$model->id} = [$model, $reason];
        } else {

            my %last_build_instdata = ( );

            my @last_build_inputs = $last_build->inputs;
            @last_build_inputs   = grep { $_->name eq 'instrument_data' } @last_build_inputs;
            %last_build_instdata = map  { $_->value_id => 1 }             @last_build_inputs;

            my @inputs = $model->instrument_data_inputs;
            my @missing_assignments_in_last_build = grep { not $last_build_instdata{$_->value_id} } @inputs;

            if (@missing_assignments_in_last_build) {
                my $reason = 'it does not have a final build with all assignments';
                $self->status_message("Requesting build of model " . $model->__display_name__ . " because $reason");
                $models_to_build{$model->id} = [$model, $reason];
            } else {
                $self->status_message("skipping rebuild of model " . $model->__display_name__ . " because all instrument data assignments are on the last build");
            }
        }
    }

    $self->status_message("Requesting builds...");

    MODEL: for my $model_and_reason (values %models_to_build) {
        my ($model, $reason) = @$model_and_reason;

        #TODO generalize via model->notify_input_build_success to make is_ready_to_build or the like
        if($model->isa('Genome::Model::SomaticValidation')) {
            if($model->tumor_sample and $model->normal_sample) {
                my @i = $model->instrument_data;

                my ($t, $n) = (0,0);
                for my $i (@i) {
                    if($i->sample eq $model->tumor_sample) { $t++; }
                    if($i->sample eq $model->normal_sample) { $n++; }
                }

                next MODEL unless ($t > 0 and $n > 0);
            }
        }

        #Will be picked up by next run of `genome model services build-queued-models`
        $model->build_requested(1, 'AQID: ' .$reason);
    }

    return 1;
}

sub _resolve_processing_for_instrument_data {
    my ($self, $instrument_data) = @_;

    my $sequencing_platform = $instrument_data->sequencing_platform;
    my $import_format = eval{ $instrument_data->import_format; };
    my @processing;
    eval {

        my $sample = $instrument_data->sample;
        unless (defined($sample)) {
            die $self->error_message('Failed to get a sample for instrument data! '.$instrument_data->id);
        }

        my $source = $sample->source; # taxon is via the source, so check it first
        if ( not $source ) {
            die $self->error_message('Failed to get a sample source for instrument data! '.$instrument_data->id);
        }

        my $taxon = $source->taxon;
        unless (defined($taxon)) {
            die $self->error_message('Failed to get a taxon from sample source for instrument data! '.$instrument_data->id);
        }

        if ($sequencing_platform eq '454') {
            if ( $instrument_data->sample->name =~ /^n\-cn?trl$/i ) { # Do not process 454 negative control (n-ctrl, n-cntrl)
                $instrument_data->ignored(1);
            }
            elsif ( $self->_is_mc16s($instrument_data) ) { # MC16s
                $self->_find_or_create_mc16s_454_qc_model($instrument_data); # always add this inst data to the QC model.
                if ( $instrument_data->read_count > 0 ) { # skip inst data w/ 0 reads
                    for my $processing_profile_id ( Genome::Model::MetagenomicComposition16s->default_processing_profile_ids ) {
                        push @processing, { processing_profile_id => $processing_profile_id, };
                    };
                }
                else {
                    $instrument_data->ignored(1);
                }
            }
            else { # skip
                $self->status_message('Skipping 454 instrument data because it is not MC16s! '.$instrument_data->id);
            }
        }
        elsif ($sequencing_platform eq 'sanger') {
            push @processing, { processing_profile_id => 2591277, };
        }
        elsif ( $import_format and $import_format eq 'genotype file' ) {
            # Genotype Microarry PP as of 2011jan25
            # Processing is not dependent on seq platform, but we use to care
            my %genotype_platforms_and_processing_ids = (
                illumina => 2166945,
                affymetrix => 2166946,
                infinium => 2575175,
                unknown => 2186707,
                plink => 2591110,
            );
            my $processing_profile_id = $genotype_platforms_and_processing_ids{$sequencing_platform};
            die $self->error_message('No genotype processing profile for platform! '.$sequencing_platform) if not $processing_profile_id;
            for my $reference_sequence_build_id (qw/ 101947881 106942997 /) {# NCBI-human-build36 => 101947881, GRCh37-lite-build37 => 106942997
                push @processing, {
                    processing_profile_id => $processing_profile_id,
                    reference_sequence_build_id => $reference_sequence_build_id,
                };
            }
        }
        elsif ($sequencing_platform eq  'solexa') {
            if($instrument_data->target_region_set_name and Genome::FeatureList->get(name => $instrument_data->target_region_set_name)->content_type eq 'validation') {
                #Do not create ref-align models--will try to assign to existing SomaticValidation models.
            } elsif ($taxon->species_latin_name =~ /homo sapiens/i) {
                if ($self->_is_pcgp($instrument_data)) {
                    push @processing, {
                        processing_profile_id => 2644306,
                        reference_sequence_build_id => 106942997,# GRCh37-lite-build37 => 106942997
                    };
                }
                elsif ( $instrument_data->sample->is_rna ) {
                    if($instrument_data->is_paired_end){
                        push @processing, {
                            processing_profile_id => $self->_default_rna_seq_processing_profile_id($instrument_data),
                            reference_sequence_build_id => 106942997,# GRCh37-lite-build37 => 106942997
                        };
                    }
                }
                else {
                    push @processing, { 
                        processing_profile_id => Genome::ProcessingProfile::ReferenceAlignment->default_profile_id,
                        reference_sequence_build_id => 106942997,# GRCh37-lite-build37 => 106942997
                    };
                }
            }
            elsif ($taxon->species_latin_name =~ /mus musculus/i){
                push @processing, {
                    processing_profile_id => Genome::ProcessingProfile::ReferenceAlignment->default_profile_id,
                    reference_sequence_build_id => 107494762,# UCSC-mouse-buildmm9 => 107494762
                };
            }
            elsif ($taxon->species_latin_name =~ /zea mays/i) {
                push @processing, {
                    processing_profile_id => Genome::ProcessingProfile::ReferenceAlignment->default_profile_id,
                    reference_sequence_build_id => 123196088,# MGSC-maize-buildB73 => 123196088
                };
            }
            elsif ($taxon->domain =~ /bacteria/i) {
                push @processing, { processing_profile_id => 2732557, };
            }
            elsif ( my @default_pp_ids = $self->_get_default_processing_profile_ids_for_instrument_data($instrument_data) ) {
                for my $processing_profile_id ( @default_pp_ids ) {
                    push @processing, { processing_profile_id => $processing_profile_id, };
                }
            }
        }

        for my $processing ( @processing ) {
            my $processing_profile_id = delete $processing->{processing_profile_id};
            my $processing_profile = Genome::ProcessingProfile->get($processing_profile_id);
            if ( not $processing_profile ) { 
                die $self->error_message('Failed to get processing profile for id! '.$processing_profile_id);
            }
            $processing->{processing_profile} = $processing_profile;

            my $reference_sequence_build_id = delete $processing->{reference_sequence_build_id};
            if ( $reference_sequence_build_id ) {
                my $reference_sequence_build = Genome::Model::Build::ReferenceSequence->get($reference_sequence_build_id);
                if ( not $reference_sequence_build ) {
                    die $self->error_message('Failed to get reference sequence build for id! '.$reference_sequence_build_id);
                }
                $processing->{reference_sequence_build} = $reference_sequence_build;
            }
        }

    };
    if($@){
        $self->error_message('Failed to get processing for instrument data id ('.$instrument_data->id.'): '.$@);
        push @processing, {error => $self->error_message};
    }
    return @processing;
}

sub _is_pcgp {
    my ($self, $instrument_data) = @_;

    foreach my $project ( $self->_get_projects_for_instrument_data($instrument_data) ) {
        my $project_id = $project->id;
        my $project_name = $project->name;

        if (defined($project_id) && grep($project_id eq $_, (2230523, 2230525, 2259255, 2342358))
            || defined($project_name) && grep($project_name =~ /$_/, ('^PCGP','^Pediatric Cancer'))) {
            return 1;
        }
    }

    return 0;
}

sub _is_mc16s {
    my ($self, $instrument_data) = @_;

    my @projects = $self->_get_projects_for_instrument_data($instrument_data);
    return if not @projects;

    my @setups = Genome::Site::TGI::Synchronize::Classes::LimsProject->get(id => [ map { $_->id } @projects ]);
    return if not @setups;

    foreach my $setup ( @setups ) {
        my $pipeline = $setup->pipeline;
        next if not $pipeline;
        return 1 if $pipeline =~ /16s/i;
    }

    return;
}

sub roi_sets {
     return (
        'WU-Space' => [
        'NCBI-human.combined-annotation-58_37c_cds_exon_and_rna_merged_by_gene',
        'NCBI-human.combined-annotation-54_36p_v2_CDSome_w_RNA',
        ],
        'TCGA-CDS' => [
        'agilent_sureselect_exome_version_2_broad_refseq_cds_only_hs37',
        'agilent sureselect exome version 2 broad refseq cds only',
        ],
    );
}

1;

