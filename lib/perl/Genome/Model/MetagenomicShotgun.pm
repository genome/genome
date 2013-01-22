package Genome::Model::MetagenomicShotgun;

use strict;
use warnings;

use Genome;
use File::Basename;
use File::Path 'make_path';
use Sys::Hostname;
use Data::Dumper;

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

    my $start_build = sub {
        my ($model, @instrument_data) = @_;
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
    };

    my $wait_build = sub {
        my $watched_build = shift;
        if ( not $watched_build ) {
            $self->status_message("No build to wait!");
            return;
        }
        $self->status_message('Watching build: '.$watched_build->__display_name__);
        $self->_wait_for_build($watched_build);
        my $status = $watched_build->status;
        if ( $status eq 'Succeeded' ) {
            $self->status_message($status.'! '.$watched_build->__display_name__);
        }
        else {
            $self->error_message($status.'! '.$watched_build->__display_name__);
        }
    };

    my $extract_data = sub {
        my ($from_build, $extraction_type) = @_;
        if ( not defined $from_build ){
            $self->status_message("No previous build provided, skipping $extraction_type data extraction");
            return;
        }
        $self->status_message("Extracting $extraction_type reads from ".$from_build->__display_name__);
        my @assignments = $from_build->instrument_data_inputs;
        my @extracted_instrument_data;
        for my $assignment (@assignments) {
            my @alignment_results = $from_build->alignment_results_for_instrument_data($assignment->value);
            if (@alignment_results > 1) {
                die $self->error_message( "multiple alignment_results found for instrument data assignment: " . $assignment->__display_name__);
            }
            if (@alignment_results == 0) {
                die $self->error_message( "no alignment_results found for instrument data assignment: " . $assignment->__display_name__);
            }
            $self->status_message("processing instrument data assignment ".$assignment->__display_name__." for unaligned reads import");

            my $alignment_result = $alignment_results[0];
            my @extracted_instrument_data_for_alignment_result = $self->_extract_data_from_alignment_result($alignment_result, $extraction_type,$self->filter_duplicates);

            push @extracted_instrument_data, \@extracted_instrument_data_for_alignment_result
        }

        unless (@extracted_instrument_data == @assignments) {
            die $self->error_message("The count of extracted instrument data sets does not match screened instrument data assignments.");
        }
        return map {@$_} @extracted_instrument_data;
    };

    my @original_instdata = $build->instrument_data;

    my $cs_build = $start_build->($contamination_screen_model, @original_instdata);
    my $cs_build_ok = $wait_build->($cs_build);
    return if not $cs_build_ok;
    my $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $cs_build, sub_model_name => 'contamination_screen');
    return if not $link_alignments;

    my @cs_unaligned;
    if ($self->filter_contaminant_fragments){
        @cs_unaligned = $extract_data->($cs_build, "unaligned paired");
    }
    else{
        @cs_unaligned = $extract_data->($cs_build, "unaligned");
    }

    my $mg_nucleotide_build = $start_build->($metagenomic_nucleotide_model, @cs_unaligned);
    my $mg_nt_build_ok = $wait_build->($mg_nucleotide_build);
    return if not $mg_nt_build_ok;
    $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_nucleotide_build, sub_model_name => 'metagenomic_nucleotide');
    return if not $link_alignments;

    my @mg_nucleotide_unaligned = $extract_data->($mg_nucleotide_build, "unaligned");
    my @mg_nucleotide_aligned = $extract_data->($mg_nucleotide_build, "aligned");

    my $mg_protein_build = $start_build->($metagenomic_protein_model, @mg_nucleotide_unaligned);
    my $mg_nr_build_ok = $wait_build->($mg_protein_build);
    return if not $mg_nr_build_ok;
    $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $mg_protein_build, sub_model_name => 'metagenomic_protein');
    return if not $link_alignments;

    my @mg_protein_aligned = $extract_data->($mg_protein_build, "aligned");

    my $viral_nucleotide_build = $start_build->($viral_nucleotide_model, @mg_nucleotide_aligned, @mg_protein_aligned);
    my $viral_protein_build = $start_build->($viral_protein_model, @mg_nucleotide_aligned, @mg_protein_aligned);

    my $viral_nr_build_ok = $wait_build->($viral_protein_build);
    return if not $viral_nr_build_ok;
    $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_protein_build, sub_model_name => 'viral_protein');
    return if not $link_alignments;

    my $viral_nt_build_ok = $wait_build->($viral_nucleotide_build);
    return if not $viral_nt_build_ok;
    $link_alignments = $self->_link_sub_build_alignments_to_build(build => $build, sub_build => $viral_nucleotide_build, sub_model_name => 'viral_nucleotide');
    return if not $link_alignments;

    return 1;
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

    $self->status_message('Find running or scheduled build for model: '.$model->__display_name__);

    UR::Context->reload('Genome::Model::Build', model_id => $model->id);

    my $build = $model->latest_build;
    if ( $build and grep { $build->status eq $_ } (qw/ Scheduled Running /) ) {
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
}

sub _extract_data_from_alignment_result{
    my ($self, $alignment, $extraction_type, $filter_duplicates) = @_;

    my $instrument_data = $alignment->instrument_data;
    my $lane = $instrument_data->lane;
    my $instrument_data_id = $instrument_data->id;

    my $aligner_name = $alignment->aligner_name; #we need to do special handling for rtg mapx alignments, which are protein space and do not produce bam files

    my $dir = $alignment->output_dir;
    my $alignment_output;
    if ($aligner_name eq 'rtg mapx'){
        unless ($extraction_type eq 'aligned'){
            die $self->error_message("metagenomic shotgun pipeline should only extract aligned reads from mapx output");
        }
        $alignment_output = $dir . '/alignments.txt';
    }else{
        $alignment_output = $dir . '/all_sequences.bam';
    }
    unless (-e $alignment_output) {
        die $self->error_message("Failed to find expected alignment output file $alignment_output\n");
    }

    # come up with the import "original_data_path", used to find existing data, and when uploading new data

    my $tmp_dir = "/tmp/extracted_reads";
    $tmp_dir .= "/".$alignment->id;
    my @instrument_data;
    my $subdir = $extraction_type;
    $subdir =~ s/\s/_/g;
    if ($filter_duplicates){
        $subdir .= "/deduplicated";
    }
    $self->status_message("Preparing imported instrument data for import path $tmp_dir/$subdir");

    my $forward_basename = "s_$lane" . "_1_sequence.txt";
    my $reverse_basename = "s_$lane" . "_2_sequence.txt";
    my $fragment_basename = "s_$lane" . "_sequence.txt";

    my $expected_data_path0 = "$tmp_dir/$subdir/$fragment_basename";
    my $expected_data_path1 = "$tmp_dir/$subdir/$forward_basename";
    my $expected_data_path2 = "$tmp_dir/$subdir/$reverse_basename";

    my $expected_se_path = $expected_data_path0;
    my $expected_pe_path = $expected_data_path1 . ',' . $expected_data_path2;

    # get any pre-existing imported instrument data for the given se_path and pe_path

    my ($se_lock, $pe_lock);

    $self->status_message("Checking for previously imported extracted reads from: $tmp_dir/$subdir");
    my $se_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_se_path);
    if ($se_instdata) {
        $self->status_message("imported instrument data already found for path $expected_se_path, skipping");
    }
    else {
        my $lock = $fragment_basename;
        $lock = $ENV{GENOME_LOCK_DIR} . '/' . $instrument_data_id . '/' . $lock;

        $self->status_message("Creating lock on $lock...");
        $se_lock = Genome::Sys->lock_resource(
            resource_lock => $lock,
            max_try => 2,
        );
        unless ($se_lock) {
            die $self->error_message("Failed to lock $expected_se_path.");
        }
    }

    my $pe_instdata;
    if ( $instrument_data->is_paired_end ){
        $pe_instdata = Genome::InstrumentData::Imported->get(original_data_path => $expected_pe_path);
        if ($pe_instdata) {
            $self->status_message("imported instrument data already found for path $expected_pe_path, skipping");
        }
        else  {
            my $lock = "s_$lane"."_1-2_sequence.txt";
            $lock = $ENV{GENOME_LOCK_DIR} . '/' . $instrument_data_id . '/' . $lock;

            $self->status_message("Creating lock on $lock...");
            $pe_lock = Genome::Sys->lock_resource(
                resource_lock => "$lock",
                max_try => 2,
            );
            unless ($pe_lock) {
                die $self->error_message("Failed to lock $expected_pe_path.");
            }
        }
    }

    if (!$se_lock and !$pe_lock) {
        $self->status_message("skipping read processing since all data is already processed and uploaded");
        return grep { defined $_ } ($se_instdata, $pe_instdata);
    }

    # extract what we did not already find...

    my $working_dir = Genome::Sys->create_temp_directory();
    my $forward_unaligned_data_path     = "$working_dir/$instrument_data_id/$forward_basename";
    my $reverse_unaligned_data_path     = "$working_dir/$instrument_data_id/$reverse_basename";
    my $fragment_unaligned_data_path    = "$working_dir/$instrument_data_id/$fragment_basename";

    my $cmd;
    if ($aligner_name eq 'rtg mapx'){
        $cmd = "genome-perl -S genome model metagenomic-shotgun mapx-alignment-to-aligned-fastq --mapx-alignment $alignment_output --output-directory $working_dir --instrument-data $instrument_data_id";
    }else{
        $cmd = "genome-perl -S gmt sam bam-to-unaligned-fastq --bam-file $alignment_output --output-directory $working_dir --ignore-bitflags"; #add ignore bitflags here because some of the aligners used in this pipeline produce untrustworthy flag information
        if ($extraction_type eq 'aligned'){
            $cmd.=" --print-aligned";
        }
        elsif($extraction_type eq 'unaligned'){
            #default
        }
        else{
            die $self->error_message("Unhandled extraction_type $extraction_type");
        }
        if ($filter_duplicates){
            $cmd.=" --filter-duplicates";
        }
    }

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        die "Failed to extract unaligned reads: $@";
    }

    my @expected_output_fastqs = ( $instrument_data->is_paired_end )
        ?  ($forward_unaligned_data_path, $reverse_unaligned_data_path, $fragment_unaligned_data_path)
        :  ($fragment_unaligned_data_path);

    my @missing = grep {! -e $_} grep { defined($_) and length($_) } @expected_output_fastqs;
    if (@missing){
        die $self->error_message(join(", ", @missing)." unaligned files missing after bam extraction");
    }
    $self->status_message("Extracted unaligned reads from bam file (@expected_output_fastqs)");

    # upload

    $self->status_message("uploading new instrument data from the post-processed unaligned reads...");
    my @properties_from_prior = qw/
        run_name
        sequencing_platform
        median_insert_size
        sd_above_insert_size
        library_name
        sample_name
    /;
    my @errors;
    my %properties_from_prior;
    for my $property_name (@properties_from_prior) {
        my $value = $instrument_data->$property_name;
        no warnings;
        $self->status_message("Value for $property_name is $value");
        $properties_from_prior{$property_name} = $value;
    }
    $properties_from_prior{subset_name} = $instrument_data->lane;

    my @source_paths = ( $instrument_data->is_paired_end )
        ?  ($fragment_unaligned_data_path, "$forward_unaligned_data_path,$reverse_unaligned_data_path" )
        :  ($fragment_unaligned_data_path);

    for my $source_data_files (@source_paths) {
        $self->status_message("Attempting to upload $source_data_files...");
        my $original_data_path;
        if ($source_data_files =~ /,/){
            $properties_from_prior{is_paired_end} = 1;
            $original_data_path = $expected_pe_path;
        }
        else {
            $properties_from_prior{is_paired_end} = 0;
            $original_data_path = $expected_se_path;
        }
        my %params = (
            %properties_from_prior,
            source_data_files => $source_data_files,
            import_format => 'sanger fastq',
        );
        $self->status_message("importing fastq with the following params:" . Data::Dumper::Dumper(\%params));

        my $command = Genome::InstrumentData::Command::Import::Fastq->create(%params);
        unless ($command) {
            $self->error_message( "Couldn't create command to import unaligned fastq instrument data!");
        };
        my $result = $command->execute();
        unless ($result) {
            die $self->error_message( "Error importing data from $source_data_files! " . Genome::InstrumentData::Command::Import::Fastq->error_message() );
        }
        $self->status_message("committing newly created imported instrument data");

        my $new_instrument_data = Genome::InstrumentData::Imported->get($command->generated_instrument_data_id);
        unless ($new_instrument_data) {
            die $self->error_message( "Failed to find new instrument data $source_data_files!");
        }

        $new_instrument_data->original_data_path($original_data_path);

        UR::Context->commit();

        if ($new_instrument_data->__changes__) {
            die "unsaved changes present on instrument data $new_instrument_data->{id} from $original_data_path!!!";
        }
        if ( $se_lock ) {
            $self->status_message("Attempting to remove lock on $se_lock...");
            unless(Genome::Sys->unlock_resource(resource_lock => $se_lock)) {
                die $self->error_message("Failed to unlock $se_lock.");
            }
            undef($se_lock);
        }
        if ( $pe_lock ) {
            $self->status_message("Attempting to remove lock on $pe_lock...");
            unless(Genome::Sys->unlock_resource(resource_lock => $pe_lock)) {
                die $self->error_message("Failed to unlock $pe_lock.");
            }
            undef($pe_lock);
        }

        push @instrument_data, $new_instrument_data;
    }

    if ( not @instrument_data ) {
        die $self->error_message("Error processing unaligned reads!");
    }

    return @instrument_data;
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
