package Genome::Site::TGI::Synchronize::UpdateApipeClasses;

use strict;
use warnings;

use Genome;
use Set::Scalar;
use Scalar::Util;
use Carp 'confess';

class Genome::Site::TGI::Synchronize::UpdateApipeClasses {
    is => 'Genome::Command::Base',
    has_optional => [
        show_object_cache_summary => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, a summary of the contents of the UR object cache is occasionally printed, useful for debugging',
        },
    ],
    has_transient_optional => [
        instrument_data_with_successful_pidfas => {
            is => 'Hash',
            default_value => {},
            doc => 'Hash of instrument data ids w/ successful PIDFAs and the path to their data file.',
        },
        _report => { is => 'Hash', default_value => {}, },
        _lock => { is => 'Text', },
    ],
    doc => 'This command contains a mapping of old LIMS-based classes to new classes that use tables in ' .
        'the MG schema and determines if anything needs to be copied over',
};

# Maps old classes to new classes. Abstract classes should not be included here because 
# it can lead to some attributes not being copied over.
sub objects_to_sync {
    return (
        'Genome::Site::TGI::Synchronize::Classes::OrganismTaxon' => 'Genome::Taxon',
        'Genome::Site::TGI::Synchronize::Classes::OrganismIndividual' => 'Genome::Individual',
        'Genome::Site::TGI::Synchronize::Classes::PopulationGroup' => 'Genome::PopulationGroup',
        'Genome::Site::TGI::Synchronize::Classes::OrganismSample' => 'Genome::Sample',
        'Genome::Site::TGI::Synchronize::Classes::LibrarySummary' => 'Genome::Library',
        'Genome::Site::TGI::Synchronize::Classes::LimsProject' => 'Genome::Project',
        'Genome::Site::TGI::Synchronize::Classes::LimsProjectSample' => 'Genome::Site::TGI::Synchronize::Classes::ProjectSample',
        'Genome::Site::TGI::Synchronize::Classes::LimsProjectInstrumentData' => 'Genome::Site::TGI::Synchronize::Classes::ProjectInstrumentData',
        'Genome::Site::TGI::Synchronize::Classes::RegionIndex454' => 'Genome::InstrumentData::454',
        'Genome::Site::TGI::Synchronize::Classes::IndexIllumina' => 'Genome::InstrumentData::Solexa',
        'Genome::Site::TGI::Synchronize::Classes::Genotyping' => 'Genome::InstrumentData::Imported',
        'Genome::Site::TGI::Synchronize::Classes::InstrumentDataAnalysisProjectBridge' => ['Genome::Config::AnalysisProject::InstrumentDataBridge', 'sync_id'],
    );
}

sub _suppress_status_messages {
    my $self = shift;

    no warnings;
    no strict 'refs';

    for my $class (qw/ 
        Genome::Model::Command::Define::Convergence
        Genome::Model::Command::Input::Update
        Genome::Model::Command::List
        Genome::ModelGroup 
        Genome::Project 
        UR::Object::Command::List
        /) {
        $class->__meta__;
        *{$class.'::status_message'} = sub{return $_[0];};
    }
    for my $class (qw/ 
        UR::Object::Command::List::Style
        /) {
        eval("use $class");
        *{$class.'::format_and_print'} = sub{return $_[0];};
    }


    return 1;
}

sub _lock_me {
    my $self = shift;
    return 1 if $ENV{UR_DBI_NO_COMMIT};
    $self->status_message('Lock...');
    my $lock = Genome::Sys->lock_resource(
        resource_lock => $ENV{GENOME_LOCK_DIR} . '/synchronize-update-apipe-classes',
        max_try => 1,
    );
    if ( not $lock ) {
        $self->error_message("Could not lock sync cron!");
        return;
    }
    $self->status_message('Lock: '.$lock);
    $self->_lock($lock);
    return 1;
}

sub _unlock_me {
    my $self = shift;
    return 1 if $ENV{UR_DBI_NO_COMMIT};
    return 1 if not $self->_lock;
    $self->status_message('Unlock...');
    eval{ Genome::Sys->unlock_resource(resource_lock => $self->_lock); };
    return 1;
}

sub execute {
    my $self = shift;
    $self->status_message('Sync LIMS to Genome...');

    $self->_lock_me;

    # Suppress overly talkative classes
    $self->_suppress_status_messages;

    # Load instrument data successful pidfas. We only sync instrument data the have a successful pidfa.
    my $load_pidfas = $self->_load_successful_pidfas;
    if ( not $load_pidfas ) {
        $self->error_message('Failed to load instruemnt data successful pidfas!');
        return;
    }

    my @classes_to_sync = $self->objects_to_sync;
    for ( my $i = 0; $i < @classes_to_sync; $i += 2 ) {
        my $lims_class = $classes_to_sync[$i];
        my $lims_ids = $self->_get_ids_for($lims_class);

        my $genome_class = $classes_to_sync[$i + 1];
        my $id_by_method;
        if (ref $genome_class eq 'ARRAY') {
            $id_by_method = $$genome_class[1];
            $genome_class = $$genome_class[0];
        }
        my $genome_ids = $self->_get_ids_for($genome_class, $id_by_method);

        $self->status_message('Detemine IDs to create...');
        my $ids_to_create = $lims_ids->difference($genome_ids);
        $self->status_message('Found IDs to create: '.scalar(@{$ids_to_create}));

        if ( not $ids_to_create->is_empty ) {
        my $start_transaction = UR::Context::Transaction->begin();

            $self->_create_genome_objects_for_lims_objects(
                ids_to_create => $ids_to_create,
                lims_class => $lims_class,
                genome_class => $genome_class,
            );
        }

        my $ids_in_genome_not_in_lims = $genome_ids->difference($lims_ids);
        $self->_report->{$genome_class}->{'missing'} = [ @$ids_in_genome_not_in_lims ];
    }

    $self->_unlock_me;

    $self->status_message('Sync LIMS to Genome...done');
    return 1;
}

sub _get_ids_for {
    my ($self, $class, $id_by_method) = @_;
    $self->status_message("Getting IDs for $class...");
    $id_by_method = 'id' unless defined($id_by_method);

    my $iterator = $class->create_iterator;
    my $set = Set::Scalar->new();
    while ( my $obj = $iterator->next ) {
        $set->insert($obj->$id_by_method);
    };
    $self->status_message("Found ".scalar(@$set)." $class IDs.");

    if ( $set->is_empty ) {
        Carp::confess('Failed to get ids for class! '.$class);
    }

    $self->status_message("Unloading $class objects...");
    my $unloaded = $class->unload;
    $self->status_message("Unloaded $unloaded objects...OK");

    $self->status_message("Getting IDs for $class...done");
    return $set;
}

sub _resolve_create_method_for {
    my ($self, $lims_class) = @_;

    Carp::confess('No LIMS class to get create method!') if not $lims_class;

    my $method_base = $lims_class;
    $method_base =~ s/Genome::Site::TGI::(Synchronize::Classes::)?//i;
    $method_base =~ s/::/_/g;
    my $create_method = '_create_' . lc($method_base);
    if ( not $self->can($create_method)) {
        Carp::confess "Did not find create method ($create_method) for LIMS class! ".$lims_class;
    }

    return $create_method;
}

sub _create_genome_objects_for_lims_objects {
    my ($self, %params) = @_;

    my $ids_to_create = delete $params{ids_to_create};
    Carp::confess('No lims ids to create genome objects!') if not $ids_to_create;
    my $lims_class = delete $params{lims_class};
    Carp::confess('No lims class to create genome objects!') if not $lims_class;
    my $genome_class = delete $params{genome_class};
    Carp::confess('No genome class to create genome objects!') if not $genome_class;

    $self->status_message('Loading LIMS objects to create in Genome...');
    my $iterator = $lims_class->create_iterator(id => [ @{$ids_to_create} ]);
    $self->status_message('Loading LIMS objects to create in Genome...done');

    my $create_method = $self->_resolve_create_method_for($lims_class);
    my $report = $self->_report;
    my $transaction = UR::Context::Transaction->begin();
    $self->status_message('Create objects in Genome...');
    while ( my $lims_obj = $iterator->next ) {
        my $genome_obj = $self->$create_method($lims_obj, $genome_class);
        if ( $genome_obj ) {
            push @{$report->{$genome_class}->{'copied'}}, $lims_obj->id;
        }
    }
    $self->status_message('Create objects in Genome...done');

    $self->status_message("Unloading $lims_class objects...");
    my $unloaded = $lims_class->unload;
    $self->status_message("Unloaded $unloaded objects...OK");

    $self->_report($report);

    $self->status_message("Commit $genome_class...");
    my $commit_ok = $transaction->commit;
    if ( not $commit_ok ) {
        $self->error_message( $transaction->error_message ) if $transaction->error_message;
        $transaction->rollback;
        Carp::confess( $self->error_message("Failed to commit $genome_class!") );

    }
    $self->status_message("Commit $genome_class...OK");

    return 1;
}

# Returns indirect and direct properties for an object and the values those properties hold
sub _get_direct_and_indirect_properties_for_object {
    my ($self, $original_object, $class, @ignore) = @_;
    my %direct_properties;
    my %indirect_properties;

    my @properties = $class->__meta__->properties;
    for my $property (@properties) {
        next if $property->is_calculated;
        next if $property->is_constant;
        next if $property->is_many;
        next if $property->id_by;
        next if $property->via and $property->via ne 'attributes';

        my $property_name = $property->property_name;
        next unless $original_object->can($property_name);
        next if @ignore and grep { $property_name eq $_ } @ignore;

        my $value = $original_object->$property_name;
        next unless defined $value;

        if ($property->via) {
            $indirect_properties{$property_name} = $value;
        }
        else {
            $direct_properties{$property_name} = $value;
        }
    }

    return (\%direct_properties, \%indirect_properties);
}

sub _load_successful_pidfas {
    my $self = shift;
    # Load successful pidfas grabbing the pidfa_output pse param, if available
    # This query/hash loading takes 10-15 secs
    print STDERR "Load instrument data successful pidfas...\n";

    my $dbh = Genome::DataSource::Oltp->get_default_handle;
    if ( not $dbh ) {
        $self->error_message('Failed to get dbh from gm schema!');
        return;
    }
    my $sql = <<SQL;
        select p1.param_value, p2.param_value
        from process_step_executions pse
        inner join pse_param p1 on p1.pse_id = pse.pse_id and p1.param_name = 'instrument_data_id'
        left join pse_param p2 on p2.pse_id = pse.pse_id and p2.param_name = 'pidfa_output'
        where pse.ps_ps_id = 3870 and pse.pr_pse_result = 'successful'
        order by p1.param_value desc
SQL

    print STDERR "PIDFA SQL:\n$sql\n";
    print STDERR "PIDFA Prepare SQL\n";
    my $sth = $dbh->prepare($sql);
    if ( not $sth ) {
        $self->error_message('Failed to prepare successful pidfa sql');
        return;
    }
    print STDERR "PIDFA Execute SQL\n";
    my $execute = $sth->execute;
    if ( not $execute ) {
        $self->error_message('Failed to execute successful pidfa sql');
        return;
    }
    print STDERR "PIDFA Fetch Results\n";
    my $instrument_data_with_successful_pidfas = $self->instrument_data_with_successful_pidfas;
    while ( my ($instrument_data_id, $pidfa_output) = $sth->fetchrow_array ) {
        # Going in reverse id order...use the most recent pidfa output for duplicate pidfas
        # pidfa output is defined for genotype microarray (genotype file) and 454 (sff file)
        $instrument_data_with_successful_pidfas->{$instrument_data_id} = $pidfa_output if not defined $instrument_data_with_successful_pidfas->{$instrument_data_id};
    }
    $sth->finish;

    print STDERR 'Loaded '.scalar(keys %$instrument_data_with_successful_pidfas)." successful PIDFAs\n";
    print STDERR 'Loaded '.scalar(grep { defined } values %$instrument_data_with_successful_pidfas)." pidfa outputs\n";
    return 1;
}

sub _create_genotyping {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required! The value is the genotype file. It must exist, too!
    my $genotype_file = $self->instrument_data_with_successful_pidfas->{$original_object->id};
    return 0 unless $genotype_file and -s $genotype_file;

    my $library_name = $original_object->sample_name.'-microarraylib';
    my ($library) = Genome::Library->get(name => $library_name, sample_id => $original_object->sample_id);
    if ( not $library ) {
        $library = Genome::Library->create(name => $library_name, sample_id => $original_object->sample_id);
        if ( not $library ) {
            Carp::confess('Failed to create genotype microarray library for sample: '.$original_object->sample_id);
        }
    }

    my $object = eval {
        $new_object_class->create(
            id => $original_object->id,
            library => $library,
            sequencing_platform => lc($original_object->platform_name),
            chip_name => $original_object->chip_name,
            version => $original_object->version,
            import_format => 'genotype file',
            import_source_name => $original_object->import_source_name,
        );
    };
    confess "Could not create new object of type $new_object_class based on object of type " .
    $original_object->class . " with id " . $original_object->id . ":\n$@" if not $object;

    my $new_genotype_file = eval{ Genome::InstrumentData::Microarray->update_genotype_file($object, $genotype_file); };
    confess "$@\nFailed to update genotype_file: $genotype_file on instrument data: ".$object->id if not $new_genotype_file;

    return 1;
}

sub _create_indexillumina {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required!
    return 0 unless exists $self->instrument_data_with_successful_pidfas->{$original_object->id};
    # Bam path required!
    return 0 unless $original_object->bam_path;

    my ($direct_properties, $indirect_properties) = $self->_get_direct_and_indirect_properties_for_object(
        $original_object,
        $new_object_class, 
        qw/ sample_name sample_id /
    );

    my $object = eval {
        $new_object_class->create(
            id => $original_object->id,
            subclass_name => $new_object_class,
            %{$direct_properties},
            %{$indirect_properties},
        );
    };
    confess "Could not create new object of type $new_object_class based on object of type " .
    $original_object->class . " with id " . $original_object->id . ":\n$@" unless $object;

    return 1;
}


sub _create_regionindex454 {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required!
    return 0 unless exists $self->instrument_data_with_successful_pidfas->{$original_object->id};

    my ($direct_properties, $indirect_properties) = $self->_get_direct_and_indirect_properties_for_object(
        $original_object,
        $new_object_class, 
        qw/ sample_name sample_id full_path/
    );
    # The value of successful pidfas is the sff file. If they are no reads, the SFF will not be defined. 
    # 454 w/ reads and no SFF should be caught in PIDFA.
    my $sff_file = $self->instrument_data_with_successful_pidfas->{$original_object->id};
    $indirect_properties->{sff_file} = $sff_file if $sff_file;

    my $object = eval {
        $new_object_class->create(
            %{$direct_properties},
            id => $original_object->id,
            subclass_name => $new_object_class,
        )
    };
    confess "Could not create new object of type $new_object_class based on object of type " .
    $original_object->class . " with id " . $original_object->id . ":\n$!" unless $object;

    my $add_attrs = $self->_add_attributes_to_instrument_data($object, $indirect_properties);
    Carp::confess('Failed to add attributes to instrument data: '.$object->__display_name__) if not $add_attrs;

    return 1;
}

sub _add_attributes_to_instrument_data {
    my ($self, $instrument_data, $attrs) = @_;

    for my $name ( keys %{$attrs} ) {
        Genome::InstrumentDataAttribute->create(
            instrument_data_id => $instrument_data->id,
            attribute_label => $name,
            attribute_value => $attrs->{$name}, 
        );
    }

    return 1;
}

sub _create_librarysummary {
    my ($self, $original_object, $new_object_class) = @_;
    return $self->_create_object($original_object, $new_object_class);
}

sub _create_organismsample {
    my ($self, $original_object, $new_object_class) = @_;
    return $self->_create_object($original_object, $new_object_class);
}

sub _create_organismindividual {
    my ($self, $original_object, $new_object_class) = @_;
    return $self->_create_object($original_object, $new_object_class);
}

sub _create_organismtaxon {
    my ($self, $original_object, $new_object_class) = @_;
    return $self->_create_object($original_object, $new_object_class);
}

sub _create_object {
    my ($self, $original_object, $new_object_class) = @_;

    my %params;
    for my $name ( $original_object->properties_to_copy ) {
        my $value = $original_object->$name;
        next if not defined $value;
        $params{$name} = $value;
    }

    my $object = eval { $new_object_class->create(%params); };
    $self->_confess_object_creation_error($original_object, $new_object_class, $@) unless $object;

    return 1;
}

sub _create_populationgroup {
    my ($self, $original_object, $new_object_class) = @_;

    # No attributes/indirect properties, etc to worry about here (except members, below)
    my %params;
    for my $property ($new_object_class->__meta__->_legacy_properties) {
        my $property_name = $property->property_name;
        $params{$property_name} = $original_object->{$property_name} if defined $original_object->{$property_name};
    }

    # Grab members from old object and pass to create parameters
    my @member_ids = $original_object->member_ids;
    if ( @member_ids ) {
        $params{member_ids} = \@member_ids;
    }

    my $object = eval {
        $new_object_class->create(
            %params, 
            id => $original_object->id, 
            subclass_name => $new_object_class
        ) 
    };
    $self->_confess_object_creation_error($original_object, $new_object_class, $@) unless $object;

    return 1;
}

sub _create_limsproject {
    my ($self, $original_object, $new_object_class) = @_;

    my $object = eval {
        $new_object_class->create(
            id => $original_object->id, 
            name => $original_object->name,
        );
    };
    $self->_confess_object_creation_error($original_object, $new_object_class, $@) unless $object;

    return 1;
}

sub _create_limsprojectinstrumentdata {
    my ($self, $original_object, $new_object_class) = @_;

    my $object = eval { 
        Genome::ProjectPart->create(
            project_id => $original_object->project_id, 
            entity_id => $original_object->instrument_data_id,
            entity_class_name => 'Genome::InstrumentData',
            label => 'instrument_data',
        );
    };
    $self->_confess_object_creation_error($original_object, $new_object_class, $@) unless $object;

    return 1;
}

sub _create_limsprojectsample {
    my ($self, $original_object, $new_object_class) = @_;

    my $object = eval { 
        Genome::ProjectPart->create(
            project_id => $original_object->project_id, 
            entity_id => $original_object->sample_id,
            entity_class_name => 'Genome::Sample',
            label => 'sample',
        );
    };
    $self->_confess_object_creation_error($original_object, $new_object_class, $@) unless $object;

    return 1;
}

sub _create_instrumentdataanalysisprojectbridge {
    my ($self, $original_object, $new_object_class) = @_;

    return $self->_create_object($original_object, $new_object_class);
}

sub _confess_object_creation_error {
    my ($self, $original_object, $new_object_class, $error) = @_;

    confess sprintf("Could not create new object of type %s based on object of type %s with id %s:\n%s",
        $new_object_class, $original_object->class, $original_object->id, $error);
}

1;

