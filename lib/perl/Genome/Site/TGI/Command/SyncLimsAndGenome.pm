package Genome::Site::TGI::Command::SyncLimsAndGenome;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();
use Genome::Utility::PluckColumn;

use Scalar::Util;
use Set::Scalar;

use constant MAX_GENOTYPE_DATA_TO_PROCESS => 500;

class Genome::Site::TGI::Command::SyncLimsAndGenome {
    is => 'Command::V2',
    has_optional => [
        expunge => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Expunge solexa and 454 instrument data from Genome that are not in LIMS.',
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
    doc => 'Copy new objects in LIMS to Genome, and remove expunged non-imported instrument data from Genome', 
};

sub execute {
    my $self = shift; 
    $self->status_message('Sync Genome from LIMS...');

    my $i_am_locked = $self->_lock_me;
    return if not $i_am_locked;

    $self->_suppress_status_messages;

    my $load_pidfas = $self->_load_successful_pidfas;
    if ( not $load_pidfas ) {
        $self->error_message('Failed to load instrument data successful pidfas!');
        return;
    }

    my $uac = $self->_update_apipe_classes;
    return if not $uac;

    if ( $self->expunge ) {
        my $expunge = $self->_expunge;
        return if not $expunge;
    }

    my $report_string = $self->_generate_report;
    print $report_string;

    $self->_unlock_me;

    $self->status_message('Sync Genome from LIMS...done');
    return 1;
}

sub _lock_me {
    my $self = shift;
    return 1 if $ENV{UR_DBI_NO_COMMIT};
    $self->status_message('Lock...');
    my $lock = Genome::Sys::LockProxy->new(
        resource => 'synchronize-genome-from-lims',
        scope => 'site',
    )->lock(
        max_try => 1,
    );
    if ( not $lock ) {
        $self->error_message("Could not lock!");
        return;
    }
    $self->status_message('Lock: ' . $lock->resource);
    $self->_lock($lock);
    return 1;
}

sub _unlock_me {
    my $self = shift;
    return 1 if $ENV{UR_DBI_NO_COMMIT};
    return 1 if not $self->_lock;
    $self->status_message('Unlock: '.$self->_lock->resource);
    eval{ $self->_lock->unlock() };
    return 1;
}

sub _suppress_status_messages {
    my $self = shift;

    no warnings;
    no strict 'refs';

    for my $class (qw/ 
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

sub _load_successful_pidfas {
    my $self = shift;
    # Load successful pidfas grabbing the pidfa_output pse param, if available
    # This query/hash loading takes 10-15 secs

    $self->status_message("Load instrument data successful pidfas...");

    my $dbh = Genome::DataSource::Oltp->get_default_handle;
    if ( not $dbh ) {
        $self->error_message('Failed to get dbh from gm schema!');
        return;
    }

    my $date_clause = '';
    unless ($self->expunge) {
        #if we don't need to check for removed data, only grab recent additions
        my $last_year = 1900 + (localtime)[5] - 1;
        $date_clause = qq{ and pse.date_scheduled > '$last_year-01-01T00:00:00Z' };
    }

    my $sql = <<SQL;
        select p1.param_value, p2.param_value
        from process_step_executions pse
        inner join pse_param p1 on p1.pse_id = pse.pse_id and p1.param_name = 'instrument_data_id'
        left join pse_param p2 on p2.pse_id = pse.pse_id and p2.param_name = 'pidfa_output'
        where pse.ps_ps_id = 3870 and pse.pr_pse_result = 'successful' $date_clause
        order by p1.param_value desc
SQL

    $self->status_message("PIDFA SQL:\n$sql");
    $self->status_message("PIDFA Prepare SQL");
    my $sth = $dbh->prepare($sql);
    if ( not $sth ) {
        $self->error_message('Failed to prepare successful pidfa sql');
        return;
    }
    $self->status_message("PIDFA Execute SQL");
    my $execute = $sth->execute;
    if ( not $execute ) {
        $self->error_message('Failed to execute successful pidfa sql');
        return;
    }
    $self->status_message("PIDFA Fetch Results");
    my $instrument_data_with_successful_pidfas = $self->instrument_data_with_successful_pidfas;
    while ( my ($instrument_data_id, $pidfa_output) = $sth->fetchrow_array ) {
        # Going in reverse id order...use the most recent pidfa output for duplicate pidfas
        # pidfa output is defined for genotype microarray (genotype file) and 454 (sff file)
        $instrument_data_with_successful_pidfas->{$instrument_data_id} = $pidfa_output if not defined $instrument_data_with_successful_pidfas->{$instrument_data_id};
    }
    $sth->finish;

    $self->status_message('Loaded '.scalar(keys %$instrument_data_with_successful_pidfas)." successful PIDFAs");
    $self->status_message('Loaded '.scalar(grep { defined } values %$instrument_data_with_successful_pidfas)." pidfa outputs");
    return 1;
}

sub _update_apipe_classes {
    my $self = shift;

    my $dictionary = Genome::Site::TGI::Synchronize::Classes::Dictionary->get;
    my @entity_names = $dictionary->entity_names;
    for my $entity_name ( @entity_names ) {
        $self->status_message("Detemine $entity_name IDs to create...");

        if ($entity_name eq 'instrument data solexa') {
            $self->_update_instrument_data($entity_name);
            next;
        }

        my $differ = Genome::Site::TGI::Command::DiffLimsAndGenome->create(
            entity_name => $entity_name,
            print_diffs => 0,
        );
        my $execute_differ = $differ->execute;
        if ( not $execute_differ ) {
            $self->error_message("Failed to execute LIMS v. Genome differ for $entity_name");
            return;
        }

        my $ids_to_create = $differ->in_lims_not_genome;
        $self->status_message('Found IDs to create: '.scalar(@{$ids_to_create}));

        my $lims_class = $differ->lims_class;
        my $genome_class_for_create = $lims_class->genome_class_for_create;
        if ( not $ids_to_create->is_empty ) {
            $self->_create_genome_objects_for_lims_objects(
                ids_to_create => $ids_to_create,
                lims_class => $lims_class,
                genome_class => $genome_class_for_create,
            );
        }

        my $in_genome_not_lims = $differ->in_genome_not_lims;
        $self->_report->{$entity_name}->{'missing'} = [ @$in_genome_not_lims ];
    }

    return 1;
}

sub _update_instrument_data {
    my ($self, $entity_name) = @_;

    my $lims_class = Genome::Site::TGI::Synchronize::Classes::Dictionary->lims_class_for_entity_name($entity_name);
    my $genome_class_for_create = $lims_class->genome_class_for_create;

    my $ids_in_genome = Genome::Utility::PluckColumn::pluck_column_from_class('Genome::InstrumentData', column_name => Genome::InstrumentData->__meta__->get_all_id_column_names);
    my @pidfa_ids = grep { Scalar::Util::looks_like_number($_) } keys %{ $self->instrument_data_with_successful_pidfas };
    my $ids_to_create = Set::Scalar->new(@pidfa_ids);
    $ids_to_create->delete(@$ids_in_genome);

    if ( not $ids_to_create->is_empty ) {
        $self->_create_genome_objects_for_lims_objects(
            ids_to_create => $ids_to_create,
            lims_class => $lims_class,
            genome_class => $genome_class_for_create,
        );
    }

    return 1;
}

sub _resolve_create_method_for {
    my ($self, $lims_class) = @_;

    Carp::confess('No LIMS class to get create method!') if not $lims_class;

    my $method_base = $lims_class;
    $method_base =~ s/Genome::Site::TGI::(Synchronize::Classes::)?//i;
    $method_base =~ s/::/_/g;
    my $create_method = '_create_' . lc($method_base);
    if ( not $self->can($create_method)) {
        $create_method = '_create_object';
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

    my $entity_name = $lims_class->entity_name;
    $self->status_message("Create $entity_name objects in Genome...");

    $self->status_message('Loading LIMS objects to create in Genome...');
    my $iterator = $lims_class->create_iterator(id => [ @{$ids_to_create} ]);
    $self->status_message('Loading LIMS objects to create in Genome...done');

    my $create_method = $self->_resolve_create_method_for($lims_class);
    my $report = $self->_report;
    my $transaction = UR::Context::Transaction->begin();
    $self->status_message('Create objects in Genome...');
    my @ids_created;
    while ( my $lims_obj = $iterator->next ) {
        my $genome_obj = $self->$create_method($lims_obj, $genome_class);
        push @ids_created, $lims_obj->id if $genome_obj;
    }
    $report->{ $lims_class->entity_name }->{copied} = \@ids_created;
    $self->status_message('Attempted: '.$ids_to_create->size);
    $self->status_message('Created:   '.@ids_created);

    $self->_report($report);

    $self->status_message("Commit $entity_name...");
    my $commit_ok = $transaction->commit;
    if ( not $commit_ok ) {
        $self->error_message( $transaction->error_message ) if $transaction->error_message;
        $transaction->rollback;
        Carp::confess( $self->error_message("Failed to commit $genome_class!") );
    }

    $self->status_message("Unloading $lims_class objects...");
    my $unloaded = $lims_class->unload;

    $self->status_message("Create $entity_name objects in Genome...done");
    return 1;
}

sub _create_object {
    my ($self, $original_object, $new_object_class) = @_;
    return $original_object->create_in_genome;
}

my $genotype_count = 0;
sub _create_genotyping {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required! The value is the genotype file. It must exist, too!
    my $genotype_file = $self->instrument_data_with_successful_pidfas->{$original_object->id};
    return 0 unless $genotype_file and -s $genotype_file;

    $genotype_count++;
    return 0 if $genotype_count > MAX_GENOTYPE_DATA_TO_PROCESS;

    $original_object->genotype_file($genotype_file);

    return $self->_create_object($original_object, $new_object_class);
}

sub _create_indexillumina {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required!
    return 0 unless exists $self->instrument_data_with_successful_pidfas->{$original_object->id};
    # Bam path required!
    return 0 unless $original_object->bam_path;

    return $self->_create_object($original_object, $new_object_class);
}

sub _create_regionindex454 {
    my ($self, $original_object, $new_object_class) = @_;

    # Successful PIDFA required!
    return 0 unless exists $self->instrument_data_with_successful_pidfas->{$original_object->id};

    my $sff_file = $self->instrument_data_with_successful_pidfas->{$original_object->id};
    $original_object->sff_file($sff_file) if $sff_file;

    return $self->_create_object($original_object, $new_object_class);
}

sub _create_instrumentdataanalysisprojectbridge {
    my ($self, $original_object, $new_object_class) = @_;

    my $inst_data = Genome::InstrumentData->get($original_object->instrument_data_id);
    return 0 unless $inst_data;
    $inst_data->unload();

    return $self->_create_object($original_object, $new_object_class);
}

sub _expunge {
    my $self = shift;

    my $report = $self->_report;

    my $dictionary = Genome::Site::TGI::Synchronize::Classes::Dictionary->get;
    for my $entity_name (keys %$report){
        my $lims_class = $dictionary->lims_class_for_entity_name($entity_name);
        my $class = $lims_class->genome_class_for_create;
        next unless $class =~ m/Genome::InstrumentData/; #only remove instrument data for now
        next if $class eq 'Genome::InstrumentData::Imported'; #imported instrument data doesn't come from LIMS, so skip it

        my @ids = @{$report->{$entity_name}->{missing}} if $report->{$entity_name}->{missing};
        next if not @ids;

        my $transaction = UR::Context::Transaction->begin();
        printf("DELETING %s %s\n", $class, join(' ', @ids));
        my @deleted;
        for my $id ( @ids ) {
            my $object = $class->get($id);
            $object->delete;
            push @deleted, $id;
        }

        my $commit_ok = $transaction->commit;
        if ( not $commit_ok ) {
            $self->error_message( $transaction->error_message ) if $transaction->error_message;
            $transaction->rollback;
            Carp::confess( $self->error_message("Failed to commit deleting missing $class!") );
        }

        $report->{$class}->{deleted} = \@deleted;
    }

    $self->_report($report);

    return 1;
}

sub _generate_report {
    my $self = shift;

    my %report = %{$self->_report};
    my $string;
    for my $type (sort keys %report) {
        $string .= "Type $type";
        for my $operation (qw/ copied missing deleted /) {
            my $num = 0;
            if (exists $report{$type}{$operation}) {
                $num = scalar @{$report{$type}{$operation}};
            }
            $string .= (', ' . (ucfirst $operation) . " $num");
        }
        $string .= "\n";
    }
    return $string;
}

1;

