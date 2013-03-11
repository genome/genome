package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;

class Genome::Model::MetagenomicShotgun {
    is => 'Genome::ModelDeprecated',
    has_param => [
        filter_contaminant_fragments => {
            is => 'Boolean',
            doc => 'when set, reads with mate mapping to contamination reference are considered contaminated, and not passed on to subsequent alignments',
            default => 0,
        },
        filter_duplicates => {
            is => 'Boolean',
            doc =>  'when set, duplicate reads are filtered out when extracting unaligned reads from contamination screen alignment',
        },
        contamination_screen_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for contamination screen',
            is_optional=> 1,
        },
        metagenomic_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for metagenomic alignment',
        },
        metagenomic_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for realignment of unaligned reads from first metagenomic alignment',
            is_optional => 1,
        },
        viral_nucleotide_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
        viral_protein_pp_id => {
            is => 'Text',
            doc => 'processing profile id to use for first viral verification alignment',
            is_optional => 1,
        },
    ],
    has => [
        contamination_screen_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'contamination_screen_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_nucleotide_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        viral_protein_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            is_optional => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'viral_protein_reference ', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        metagenomic_nucleotide_reference => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_mutable => 1,
            via => 'inputs',
            to => 'value',
            where => [name => 'metagenomic_alignment_reference', value_class_name => 'Genome::Model::Build::ImportedReferenceSequence'],
        },
        contamination_screen_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'contamination_screen_model'],
        },
        metagenomic_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            via => 'from_model_links', 
            to => 'from_model',
            where => [role => 'metagenomic_nucleotide_model'],
        },
        metagenomic_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1,
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'metagenomic_protein_model'],
        },
        viral_nucleotide_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_nucleotide_model'],
        },
        viral_protein_model => {
            is => 'Genome::Model::ReferenceAlignment',
            is_optional => 1, 
            via => 'from_model_links',
            to => 'from_model',
            where => [role => 'viral_protein_model'],
        },
    ],
};

sub sub_model_labels {
    return (qw/ contamination_screen metagenomic_nucleotide metagenomic_protein viral_nucleotide viral_protein /);
}

sub build_subclass_name {
    return 'metagenomic-composition-shotgun';
}

sub delete {
    my $self = shift;
    for my $sub_model ($self->from_models) {
        $sub_model->delete;
    }
    return $self->SUPER::delete(@_);
}

sub create {
    my ($class, %params) = @_;

    my $self = $class->SUPER::create(%params);
    return unless $self;

    my $processing_profile = $self->processing_profile;
    for ( $self->sub_model_labels ){
        my $pp_method = $_."_pp_id";
        if($self->processing_profile->$pp_method) {
            my $model = $self->_create_model_for_type($_);
            unless ($model) {
                $self->error_message("Error creating $_ model!");
                $self->delete;
                return;
            }
        }
    }

    return $self;
}

sub sequencing_platform{
    return 'solexa';
}

# TODO make this get or create
sub _create_model_for_type {
    my $self = shift;
    my $type = shift;

    #CREATE UNDERLYING REFERENCE ALIGNMENT MODELS
    my $pp_accessor = $type."_pp_id";
    my $reference_accessor = $type."_reference";
    my %model_params = (
        processing_profile_id => $self->processing_profile->$pp_accessor,
        subject_name => $self->subject_name,
        name => $self->name.".$type model",
        reference_sequence_build=>$self->$reference_accessor,
    );
    my $model = Genome::Model::ReferenceAlignment->create( %model_params );

    unless ($model){
        die $self->error_message("Couldn't create contamination screen model with params ".join(", ", map {$_ ."=>". $model_params{$_}} keys %model_params) );
    }

    $self->add_from_model(from_model=> $model, role=>$type.'_model');
    $self->status_message("Created $type model ".$model->__display_name__);

    return $model;
}

sub _resolve_resource_requirements_for_build {
    my ($self, $build) = @_;
    my @instrument_data = $build->instrument_data;
    my $gtmp = 30 + 5 * (1 + scalar(@instrument_data));
    return "-R 'rusage[gtmp=$gtmp:mem=16000]' -M 16000000";
}

sub _execute_build {
    my ($self, $build) = @_;

    my $model = $build->model;
    $self->status_message('Build '.$model->__display_name__);
    my $contamination_screen_model = $model->contamination_screen_model;
    $self->status_message("Got contamination_screen_model ".$contamination_screen_model->__display_name__) if $contamination_screen_model;
    my $metagenomic_nucleotide_model = $model->metagenomic_nucleotide_model;
    $self->status_message("Got metagenomic_nucleotide_model ".$metagenomic_nucleotide_model->__display_name__) if $metagenomic_nucleotide_model;
    my $metagenomic_protein_model = $model->metagenomic_protein_model;
    $self->status_message("Got metagenomic_protein_model ".$metagenomic_protein_model->__display_name__) if $metagenomic_protein_model;
    my $viral_nucleotide_model = $model->viral_nucleotide_model;
    $self->status_message("Got viral_nucleotide_model ".$viral_nucleotide_model->__display_name__) if $viral_nucleotide_model;
    my $viral_protein_model = $model->viral_protein_model;
    $self->status_message("Got viral_protein_model ".$viral_protein_model->__display_name__) if $viral_protein_model;

    my @cs_unaligned = $self->_screen_contaminants($build);
    return if not @cs_unaligned;

    my %mg_nucleotide_results = $self->_run_meta_nt($build, @cs_unaligned);
    return if not %mg_nucleotide_results;

    my @mg_nucleotide_unaligned = @{$mg_nucleotide_results{unaligned}};
    my @mg_protein_aligned = $self->_run_meta_nr($build,  @mg_nucleotide_unaligned);
    return if not @mg_protein_aligned;

    my @mg_nucleotide_aligned = @{$mg_nucleotide_results{aligned}};
    my $viral_protein_build = $self->_start_build($viral_protein_model, @mg_nucleotide_aligned, @mg_protein_aligned);
    my $viral_nr_build_ok = $self->_wait_for_build($viral_protein_build);
    return if not $viral_nr_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_protein_build, sub_model_name => 'viral_protein');
    return if not $link_alignments;

    my $viral_nucleotide_build = $self->_start_build($viral_nucleotide_model, @mg_nucleotide_aligned, @mg_protein_aligned);
    my $viral_nt_build_ok = $self->_wait_for_build($viral_nucleotide_build);
    return if not $viral_nt_build_ok;
    $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_nucleotide_build, sub_model_name => 'viral_nucleotide');
    return if not $link_alignments;

    return 1;
}

sub _screen_contaminants {
    my ($self, $build) = @_;

    my $contamination_screen_model = $self->contamination_screen_model;
    my @original_instdata = $build->instrument_data;
    my $cs_build = $self->_start_build($contamination_screen_model, @original_instdata);
    my $cs_build_ok = $self->_wait_for_build($cs_build);
    return if not $cs_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $cs_build, sub_model_name => 'contamination_screen');
    return if not $link_alignments;

    my @cs_unaligned;
    if ($self->filter_contaminant_fragments){
        @cs_unaligned = $self->_extract_data($cs_build, "unaligned paired");
    }
    else{
        @cs_unaligned = $self->_extract_data($cs_build, "unaligned");
    }

    return @cs_unaligned;
}

sub _run_meta_nt {
    my ($self, $build, @cs_unaligned) = @_;

    my $metagenomic_nucleotide_model = $build->model->metagenomic_nucleotide_model;
    my $mg_nucleotide_build = $self->_start_build($metagenomic_nucleotide_model, @cs_unaligned);
    my $mg_nt_build_ok = $self->_wait_for_build($mg_nucleotide_build);
    return if not $mg_nt_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_nucleotide_build, sub_model_name => 'metagenomic_nucleotide');
    return if not $link_alignments;

    my @mg_nucleotide_aligned = $self->extract_data($mg_nucleotide_build, "aligned");
    my @mg_nucleotide_unaligned = $self->_extract_data($mg_nucleotide_build, "unaligned");

    return ( 
        aligned => \@mg_nucleotide_aligned,
        unaligned => \@mg_nucleotide_unaligned
    );
}

sub _run_meta_nr {
    my ($self, $build, @mg_nucleotide_unaligned) = @_;

    my $metagenomic_protein_model = $build->model->metagenomic_protein_model;
    my $mg_protein_build = $self->_start_build($metagenomic_protein_model, @mg_nucleotide_unaligned);
    my $mg_nr_build_ok = $self->_wait_for_build($mg_protein_build);
    return if not $mg_nr_build_ok;

    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_protein_build, sub_model_name => 'metagenomic_protein');
    return if not $link_alignments;

    my @mg_protein_aligned = $self->_extract_data->($mg_protein_build, "aligned");

    return @mg_protein_aligned;
}

sub _start_build  {
    my ($self, $model, @instrument_data) = @_;

    my %existing;
    for my $inst_data($model->instrument_data){
        $existing{$inst_data->id} = $inst_data;
    }
    my @to_add;
    for my $inst_data(@instrument_data){
        if ($existing{$inst_data->id}) {
            delete $existing{$inst_data->id};
        }
        else {
            push @to_add, $inst_data;
        }
    }
    for my $inst_data (values %existing){
        $self->status_message("Removing Instrument Data " . $inst_data->id . " from model " . $model->__display_name__);
        $model->remove_instrument_data($inst_data)
    }
    for my $inst_data(@to_add){
        $self->status_message("Adding Instrument Data " . $inst_data->id . " to model " . $model->__display_name__);
        $model->add_instrument_data($inst_data);
    }
    if (@instrument_data == 0){
        $self->status_message("No instrument data for model ".$model->__display_name__.", skipping build");
        return;
    }

    my $build = $self->_build_if_necessary($model);
    return $build;
}

sub _build_if_necessary {
    my ($self, @models) = @_;

    my (@succeeded_builds, @watched_builds);
    for my $model ( @models ) {
        $self->status_message('Model: '. $model->__display_name__);
        $self->status_message('Search for succeeded build');
        my $succeeded_build = $model->last_succeeded_build;
        if ( $succeeded_build and $self->_verify_model_and_build_instrument_data_match($model, $succeeded_build) ) {
            $self->status_message('Found succeeded build: '.$succeeded_build->__display_name__);
            push @succeeded_builds, $succeeded_build;
            next;
        }
        $self->status_message('No succeeded build');
        $self->status_message('Search for scheduled or running build');
        my $watched_build = $self->_find_scheduled_or_running_build_for_model($model);
        if ( not $watched_build ) {
            $self->status_message('No scheduled or running build');
            $self->status_message('Start build');
            $watched_build = $self->_start_build_for_model($model);
            return if not $watched_build;
        }
        $self->status_message('Watching build: '.$watched_build->__display_name__);
        push @watched_builds, $watched_build;
    }

    my @builds = (@succeeded_builds, @watched_builds);
    if ( not @builds ) {
        $self->error_message('Failed to find or start any builds');
        return;
    }

    if ( @models != @builds ) {
        $self->error_message('Failed to find or start a build for each model');
        return;
    }

    return ( @builds > 1 ? @builds : $builds[0] );
}

sub _verify_model_and_build_instrument_data_match {
    my ($self, $model, $build) = @_;

    Carp::confess('No model to verify instrument data') if not $model;
    Carp::confess('No build to verify instrument data') if not $build;

    my @build_instrument_data = sort {$a->id <=> $b->id} $build->instrument_data;
    my @model_instrument_data = sort {$a->id <=> $b->id} $model->instrument_data;

    $self->status_message('Model: '.$model->__display_name__);
    $self->status_message('Model instrument data: '.join(' ', map { $_->id } @model_instrument_data));
    $self->status_message('Build: '.$build->__display_name__);
    $self->status_message('Build instrument data: '.join(' ', map { $_->id } @build_instrument_data));

    if ( @build_instrument_data != @model_instrument_data ) {
        $self->status_message('Model and build instrument data count does not match');
        return;
    }

    for ( my $i = 0; $i < @model_instrument_data; $i++ ) {
        my $build_instrument_data = $build_instrument_data[$i];
        my $model_instrument_data = $model_instrument_data[$i];

        if ($build_instrument_data->id ne $model_instrument_data->id) {
            $self->status_message("Missing instrument data.");
            return;
        }
    }

    return 1;
}

sub _find_scheduled_or_running_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('No model sent to find running or scheduled build') if not $model;

    $self->status_message('Looking for running or scheduled build for model: '.$model->__display_name__);

    UR::Context->reload('Genome::Model', id => $model->id);
    UR::Context->reload('Genome::Model::Build', model_id => $model->id);
    UR::Context->reload('Genome::Model::Event', model_id => $model->id);

    my $build = $model->latest_build;
    return if not $build;

    $self->status_message( sprintf('Build: %s %s', $build->id, $build->status) );
    if ( grep { $build->status eq $_ } (qw/ Scheduled Running /) ) {
        return $build;
    }

    return;
}

sub _start_build_for_model {
    my ($self, $model) = @_;

    Carp::confess('no model sent to start build') if not $model;

    my $cmd = 'genome model build start '.$model->id.' --job-dispatch apipe --server-dispatch workflow'; # these are defaults
    $self->status_message('cmd: '.$cmd);

    UR::Context->commit();
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        die $self->error_message('failed to execute build start command');
    }

    my $build = $self->_find_scheduled_or_running_build_for_model($model);
    if ( not $build ) {
        die $self->error_message('executed build start command, but cannot find build.');
    }

    return $build;
}

sub _wait_for_build {
    my ($self, $build) = @_;

    if ( not $build ) {
        $self->status_message("No build to wait!");
        return;
    }
    $self->status_message('Watching build: '.$build->__display_name__);

    my $last_status = '';
    my $time = 0;
    my $inc = 30;
    while (1) {
        UR::Context->current->reload($build->the_master_event);
        my $status = $build->status;
        if ($status and !($status eq 'Running' or $status eq 'Scheduled')){
            return 1;
        }

        if ($last_status ne $status or !($time % 300)){
            $self->status_message("Waiting for build(~$time sec) ".$build->id.", status: $status");
        }
        sleep $inc;
        $time += $inc;
        $last_status = $status;
    }

    my $status = $build->status;
    if ( $status eq 'Succeeded' ) {
        $self->status_message($status.'! '.$build->__display_name__);
        return 1;
    }
    else {
        $self->error_message($status.'! '.$build->__display_name__);
        return;
    }
}

sub _extract_data {
    my ($self, $from_build, $extraction_type) = @_;


    my $extract_from_alignment = Genome::Model::MetagenomicShotgun::Build::ExtractFromAlignment->create(
        build => $self,
        sub_build => $from_build,
        type => $extraction_type,
    );
    if ( not $extract_from_alignment ) {
        $self->error_message("Failed to create extract from alignment! Tried to extract $extraction_type from ".$from_build->__display_name__);
        return;
    }

    my $execute_ok = $extract_from_alignment->execute;
    if ( not $execute_ok ) {
        $self->error_message("Failed to execute extract from alignment! Tried to extract $extraction_type from ".$from_build->__display_name__);
        return;
    }

    return 1;
}

sub _link_sub_build_alignments_to_build {
    my ($self, %params) = @_;

    my $build = delete $params{build};
    Carp::confess('No build given to link alignments!') if not $build;
    my $sub_build = delete $params{sub_build};
    Carp::confess('No sub-build given to link alignments!') if not $sub_build;
    my $sub_model_name = delete $params{sub_model_name};
    Carp::confess('No sub model name given to link alignments!') if not $sub_model_name;
    Carp::confess('Unknown params given to _link_sub_build_alignments_to_build! '.Data::Dumper::Dumper(\%params)) if %params;

    my $dir = $build->data_directory;
    my $sub_dir = $dir.'/'.$sub_model_name;
    my $create_ok = eval{ Genome::Sys->create_directory($sub_dir); };
    if ( not $create_ok ) {
        $self->error_message($@) if $@;
        $self->error_message("Failed to create $sub_model_name sub dir! ".$sub_dir);
        return;
    }

    for my $instrument_data ( $sub_build->instrument_data ) {
        my @alignments = $sub_build->alignment_results_for_instrument_data($instrument_data); # This should only be one.
        for my $alignment ( @alignments ) {
            my $target = $alignment->output_dir;
            my $link = $sub_dir.'/'.$instrument_data->id.'-'.$alignment->id;
            unlink $link;
            my $link_ok = eval{ Genome::Sys->create_symlink($target, $link); };
            if ( not $link_ok ) {
                $self->error_message($@) if $@;
                $self->error_message("Failed to create symlink! From $link to $target");
                return;
            }
        }
    }

    return 1;
}

1;
