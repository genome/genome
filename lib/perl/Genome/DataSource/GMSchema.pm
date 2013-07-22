package Genome::DataSource::GMSchema; # an exact copy of ::Main except the name

use strict;
use warnings;
use Genome;
use Carp;
use File::lockf;
use DBD::Pg;
use IO::Socket;

class Genome::DataSource::GMSchema {
    is => 'UR::DataSource::Pg',
    has_constant => [
        server => { default_value => 'dbname=genome' },
        login => { default_value => 'genome' },
        auth => { default_value => 'changeme' },
        owner => { default_value => 'public' },
    ],
};

sub _ds_tag {
    'Genome::DataSource::GMSchema';
}

sub clone_db_handles_for_child_process {
    my $self = shift;

    Genome::DataSource::GMSchemaOracle->get()->clone_db_handles_for_child_process;
    Genome::DataSource::PGTest->get()->clone_db_handles_for_child_process;

    return $self->SUPER::clone_db_handles_for_child_process;
}

our $THIS_COMMIT_ID = 'not within _sync_database';

# This datasource now commits to both Oracle AND postgres. The postgres commit is
# done within an eval so its result does not in any way affect the Oracle commit.
sub _sync_database {
    my $self = shift;
    my %params = @_;

    # Need to remove all commit observers, they will be fired during commit and that's no good!
    #my @observers = UR::Observer->get();
    #for my $observer (@observers) {
    #    $observer->delete;
    #}

    # Need to update all classes that have been changed (and all of their parent
    # classes) to use the postgres datasource instead of Oracle.
    my %classes = map { $_->class => 1 } @{$params{changed_objects}};
    for my $class (sort keys %classes) {
        my $meta = UR::Object::Type->get($class);
        next unless $meta;
        my @metas = ($meta, $meta->ancestry_class_metas);
        for my $meta (@metas) {

            # If the object has been deleted and we're dealing with the Ghost, need to find the non-Ghost class
            # and add that to the list of metas we need to update. Otherwise, that class won't get updated to use
            # postgres and an error will occur.
            if ($meta->class_name =~ /::Ghost$/) {
                my $non_ghost_class = $meta->class_name;
                $non_ghost_class =~ s/::Ghost$//;
                my $non_ghost_meta = UR::Object::Type->get($non_ghost_class);
                if ($non_ghost_meta) {
                    push @metas, $non_ghost_meta;
                }
                else {
                    Carp::confess "Could not find meta object for non-ghost class $non_ghost_class!";
                }
            }

            if (defined $meta->data_source) {
                $self->rewrite_classdef_to_use_postgres($meta);     
            }
        }
    }

    # Update meta datasource to point to an empty file so we don't get failures due to not being
    # able to find column/table meta data.
    my $temp_file_fh = File::Temp->new;
    my $temp_file = $temp_file_fh->filename;
    $temp_file_fh->close;

    my $meta_ds = Genome::DataSource::Meta->_singleton_object;
    $meta_ds->get_default_handle->disconnect;
    $meta_ds->server($temp_file);
    
    return $self->SUPER::_sync_database(@_);
}

sub log_error {
    my $error = shift;
    my $log_string = create_log_message($error);
    my $log_fh = open_error_log();
    # If we can't get a file handle to the log file, no worries. Just continue without making a peep.
    if ($log_fh) {
        my $lock_status = File::lockf::lock($log_fh);
        # this returns 0 on success
        unless ($lock_status != 0) {
            $log_fh->print("$log_string\n");
            File::lockf::ulock($log_fh);
        }
        $log_fh->close;
    }
}

sub log_commit_time {
    my($db_name, $time) = @_;

    my $original_cmdline = get_command_line();
    my $execution_id = $ENV{'GENOME_EXECUTION_ID'} || '';

    my $path = _determine_base_log_pathname();
    $path .= "-${db_name}.timing";
    my $fh = IO::File->new($path,'>>');
    unless ($fh) {
        print STDERR "Couldn't open $path for appending: $!\n";
        return;
    }
    chmod(0666,$path);
    my $lock_status = File::lockf::lock($fh);
    # lock gives us a zero return value on success
    unless ($lock_status != 0) {
        $fh->print(join("\t",$THIS_COMMIT_ID, $execution_id, $time, $original_cmdline) . "\n");
        File::lockf::ulock($fh);
    }
    $fh->close();
}

# called whenever we generate an ID
sub autogenerate_new_object_id_for_class_name_and_rule {
    my $self = shift;
    UR::Object::Type->autogenerate_new_object_id_uuid; 
}

# called before attempting to write an SQL query
our %rewritten;
sub _generate_class_data_for_loading {
    my ($self, $meta) = @_;
    #$DB::single = 1 if $meta->{class_name} eq 'Genome::SubjectAttribute';
    my @ancestor_metas = $meta->ancestry_class_metas;
    for my $m ($meta, @ancestor_metas) {
        my $ds = $m->data_source;
        if ($ds and $ds == $self) {
            $self->rewrite_classdef_to_use_postgres($m);
        }
    }
    $self->SUPER::_generate_class_data_for_loading($meta);
}

# used before query _and_ from _sync_database
sub rewrite_classdef_to_use_postgres {
    my $self = shift;
    my $meta = shift;

    if ($rewritten{$meta->id}) {
        return;
    }
    $rewritten{$meta->id} = 1;

    if (my ($entity_class) = ($meta->id =~ /^(.*)::Ghost$/)) {
        my $other_meta = $entity_class->__meta__;
        return $self->rewrite_classdef_to_use_postgres($other_meta);
    }
    
    my $class = $meta->class_name;

    $meta->data_source_id($self->_ds_tag);

    # Columns are stored directly on the meta object as an optimization, need to be updated
    # in addition to table/column objects.
    my (undef, $cols) = @{$meta->{'_all_properties_columns'}};
    $_ = lc $_ foreach (@$cols);

    if (defined $meta->table_name) {
        my $oracle_table = $meta->table_name;
        my $postgres_table = $self->postgres_table_name_for_oracle_table($oracle_table);
        unless ($postgres_table) {
            Carp::confess "Could not find postgres equivalent for oracle table $oracle_table while working on class $class";
        }
        $meta->table_name($postgres_table);

        my @properties = $meta->all_property_metas;
        for my $property (@properties) {
            next unless $property->column_name;
            $property->column_name(lc $property->column_name);
        }
    }
}

sub postgres_table_name_for_oracle_table {
    my $self = shift;
    my $oracle_table = shift;
    return unless $oracle_table;
    my %mapping = $self->oracle_to_postgres_table_mapping;
    return $mapping{lc $oracle_table};
}

sub OLD_pg {
    my $self = shift;
    return 1 unless $ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE};

    my @o = grep { ref($_) eq 'UR::DeletedRef' } UR::Context->all_objects_loaded('UR::Object');
    if (@o) {
        print Data::Dumper::Dumper(\@o);
        Carp::confess();
    }

    # Determine what has changed.
    my @changed_objects = (
        UR::Context->all_objects_loaded('UR::Object::Ghost'),
        grep { $_->__changes__ } UR::Context->all_objects_loaded('UR::Object')
        #UR::Util->mapreduce_grep(sub { $_[0]->__changes__ },$self->all_objects_loaded('UR::Object'))
    );

    my @real_changed_objects = grep {UR::Context->resolve_data_source_for_object($_)} @changed_objects;

    return 1 unless (@real_changed_objects);


    print "Database updating has been paused, please wait until updating has been resumed...\n";

    my @data_sources = UR::Context->all_objects_loaded('UR::DataSource::RDBMS');
    for my $ds (@data_sources) {
        $ds->disconnect_default_handle if $ds->has_default_handle;
    }

    while (1) {
        sleep sleep_length();
        last unless $ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE};
    }

    print "Database updating has been resumed, continuing commit!\n";
    return 1;
}

sub oracle_to_postgres_table_mapping {
    return (
        'search_index_queue' => 'web.search_index_queue',
        'feature_list' => 'model.feature_list',
        'fragment_library' => 'instrument.fragment_library',
        'genome_disk_allocation' => 'disk.allocation',
        'disk_volume' => 'disk.volume',
        'disk_group' => 'disk.group',
        'disk_volume_group' => 'disk.volume_group_bridge',
        'genome_model' => 'model.model',
        'genome_model_build' => 'model.build',
        'genome_model_build_input' => 'model.build_input',
        'genome_model_build_link' => 'model.build_link',
        'genome_model_metric' => 'model.build_metric',
        'genome_model_event' => 'model.event',
        'genome_model_event_input' => 'model.event_input',
        'genome_model_event_metric' => 'model.event_metric',
        'genome_model_event_output' => 'model.event_output',
        'genome_model_group' => 'model.model_group_bridge',
        'genome_model_input' => 'model.model_input',
        'genome_model_link' => 'model.model_link',
        'genome_nomenclature' => 'web.nomenclature',
        'genome_nomenclature_enum_value' => 'web.nomenclature_enum_value',
        'genome_nomenclature_field' => 'web.nomenclature_field',
        'genome_project' => 'subject.project',
        'genome_project_part' => 'subject.project_part',
        'genome_subject' => 'subject.subject',
        'genome_subject_attribute' => 'subject.subject_attribute',
        'genome_sys_user' => 'subject.user',
        'genome_sys_user_role' => 'subject.role',
        'genome_task' => 'web.task',
        'genome_task_params' => 'web.task_params',
        'instrument_data' => 'instrument.data',
        'instrument_data_attribute' => 'instrument.data_attribute',
        'misc_attribute' => 'subject.misc_attribute',
        'misc_note' => 'subject.misc_note',
        'model_group' => 'model.model_group',
        'processing_profile' => 'model.processing_profile',
        'processing_profile_param' => 'model.processing_profile_param',
        'software_result' => 'result.software_result',
        'software_result_input' => 'result.input',
        'software_result_metric' => 'result.metric',
        'software_result_param' => 'result.param',
        'software_result_user' => 'result.user',
        'genome_model_build_variant' => 'model.build_variant',
        'genome_model_variant' => 'model.variant',
    );
}

1;
__END__
sub sleep_length {
    return 30;
}

sub create_iterator_closure_for_rule {
  my $self = shift;

  $self->pause_db_queries_if_necessary();

  $self->SUPER::create_iterator_closure_for_rule(@_);
}

sub create_dbh {
  my $self = shift;
  $self->pause_db_queries_if_necessary();


  $self->SUPER::create_dbh(@_);
}

sub pause_db_queries_if_necessary {
    my $self = shift;
    return 1 unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};

    print "Database querying has been paused; disconnecting db handles.  Please wait until the query pause is released...\n";

    my @data_sources = UR::Context->all_objects_loaded('UR::DataSource::RDBMS');
    for my $ds (@data_sources) {
        $ds->disconnect_default_handle if $ds->has_default_handle;
    }

    while (1) {
        sleep 30;
        last unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};
    }

    print "Database updating has been resumed, continuing query!\n";
    return 1;
}

# A list of the old GM schema tables that Genome::Model should ignore
my @OLD_GM_TABLES = qw(
ALL_ALLELE_TYPE
ANALYSIS_METHOD
CHROMOSOME
COLLABORATOR
COLLABORATOR_SAMPLE
CONSERVATION_SCORE
EGI_TYPE
EXTERNAL_GENE_ID
GENE
GENE_EXPRESSION
GENE_GENE_EXPRESSION
GENE_GENOMIC_FEATURE
GENE_GENOTYPE
GENOME_MODEL_1
GENOME_UPDATE_HISTORY
GENOMIC_FEATURE
GENOTYPE
GENOTYPE_VARIATION
GE_DETECTION
GE_TECH_TYPE
GF_FEATURE_TYPE
GO_XREF
GROUP_INFO
GROUP_TYPES
HISTOLOGY
IPRO_GENE_TRANSCRIPT_XREF
IPRO_RESULTS
MAF
META_GROUP
META_GROUP_STATUS
PF_FEATURE_TYPE
PP_DETECTION_SOFT
PP_MAPPING_REFERENCE
PP_TECH_TYPE
PROCESS_PROFILE
PROCESS_PROFILE_SOURCE
PROTEIN
PROTEIN_FEATURE
PROTEIN_VARIATION
PROTEIN_VARIATION_SCORE
PVS_SCORE_TYPE
PVS_SOFTWARE
READ_GROUP
READ_GROUP_GENOTYPE
READ_GROUP_GENOTYPE_OLD
READ_GROUP_INFO
READ_INFO
REPEAT_INFO
RGGI_INFO_TYPE
RGG_INFO
RGG_INFO_OLD
RI_REPEAT_TYPE
SAMPLE
SAMPLE_GENE
SAMPLE_GENE_EXPRESSION
SAMPLE_GENOTYPE
SAMPLE_GROUP_INFO
SAMPLE_HISTOLOGY
SAMPLE_SUBTYPE
SAMPLE_TISSUE
SAMPLE_TYPE
SEQUENCE_DIFF
SEQUENCE_DIFF_EVAL
SEQUENCE_DIFF_PART
SGIAM_SCORE_TYPE
SGI_ANALYSIS_METHOD
SG_INFO_TYPE
STRAND
SUBMITTER
SUBMITTER_METHOD
TERM
TISSUE
TRANSCRIPT
TRANSCRIPT_SOURCE
TRANSCRIPT_STATUS
TRANSCRIPT_SUB_STRUCTURE
TRANSCRIPT_VARIATION
TSS_STRUCTURE_TYPE
TV_OLD
TV_TYPE
VARIANT_REVIEW_DETAIL_OLD
VARIANT_REVIEW_LIST_FILTER
VARIANT_REVIEW_LIST_MEMBER
VARIATION
VARIATION_FREQUENCY
VARIATION_GROUP
VARIATION_INSTANCE
VARIATION_ORIG
VARIATION_ORIG_VARIATION_TYPE
VARIATION_SCORE
VS_SCORE_TYPE
VS_SOFTWARE
);

sub _ignore_table {
    my($self,$table_name) = @_;

    return 1 if $self->SUPER::_ignore_table($table_name);

    return scalar(grep { $_ eq $table_name } @OLD_GM_TABLES);
}

1;

